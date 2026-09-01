import argparse
import json
import re
from pathlib import Path

import h5py
import numpy as np
import pandas as pd


def decode(values):
    return np.array(
        [
            value.decode("utf-8", "replace")
            if isinstance(value, (bytes, np.bytes_))
            else str(value)
            for value in values
        ],
        dtype=object,
    )


def categorical(group):
    return decode(group["categories"][:]), group["codes"][:].astype(np.int32)


def slugify(value):
    return re.sub(r"[^a-z0-9]+", "_", value.lower()).strip("_")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--h5ad", required=True)
    parser.add_argument("--pseudobulk-metadata", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--min-cells", type=int, default=10)
    parser.add_argument("--gene-block-size", type=int, default=64)
    args = parser.parse_args()

    h5ad_path = Path(args.h5ad)
    out_dir = Path(args.out)
    count_dir = out_dir / "count_matrices"
    metadata_dir = out_dir / "sample_metadata"
    for directory in (out_dir, count_dir, metadata_dir):
        directory.mkdir(parents=True, exist_ok=True)

    metadata = pd.read_csv(
        args.pseudobulk_metadata,
        dtype={"donor_id": str, "sample_id": str},
    )
    metadata = metadata[metadata["cell_count"] >= args.min_cells].copy()
    metadata = metadata.sort_values(
        ["cell_type", "donor_id", "seq_batch"], kind="stable"
    ).reset_index(drop=True)
    if metadata.duplicated(["cell_type", "sample_id"]).any():
        raise ValueError("Duplicate cell-type and sample_id combinations")
    metadata["group_index"] = np.arange(len(metadata), dtype=np.int32)

    with h5py.File(h5ad_path, "r") as handle:
        donor_names, donor_codes = categorical(handle["obs/donor_id"])
        cell_type_names, cell_type_codes = categorical(handle["obs/cell_type"])
        seq_batch_names, seq_batch_codes = categorical(handle["obs/seq_batch"])
        donor_lookup = {value: index for index, value in enumerate(donor_names)}
        cell_type_lookup = {
            value: index for index, value in enumerate(cell_type_names)
        }
        seq_batch_lookup = {
            value: index for index, value in enumerate(seq_batch_names)
        }

        n_full_groups = len(donor_names) * len(cell_type_names) * len(seq_batch_names)
        full_to_compact = np.full(n_full_groups, -1, dtype=np.int32)
        for row in metadata.itertuples(index=False):
            full_code = (
                donor_lookup[row.donor_id]
                * (len(cell_type_names) * len(seq_batch_names))
                + cell_type_lookup[row.cell_type] * len(seq_batch_names)
                + seq_batch_lookup[row.seq_batch]
            )
            full_to_compact[full_code] = row.group_index

        cell_full_groups = (
            donor_codes * (len(cell_type_names) * len(seq_batch_names))
            + cell_type_codes * len(seq_batch_names)
            + seq_batch_codes
        )
        cell_to_group = full_to_compact[cell_full_groups]

        gene_ids = decode(handle["raw/var/_index"][:])
        gene_symbols = decode(handle["raw/var/name"][:])
        if len(set(gene_symbols)) != len(gene_symbols):
            raise ValueError("Duplicate gene symbols in raw/var/name")
        genes = pd.DataFrame(
            {
                "gene_index": np.arange(len(gene_ids), dtype=np.int32),
                "gene_id": gene_ids,
                "gene_symbol": gene_symbols,
            }
        )
        genes.to_csv(out_dir / "genes.csv", index=False)

        n_genes = len(gene_ids)
        n_groups = len(metadata)
        all_counts_path = out_dir / "all_groups_gene_by_group.int32.bin"
        all_counts = np.memmap(
            all_counts_path,
            dtype=np.int32,
            mode="w+",
            shape=(n_genes, n_groups),
            order="C",
        )

        raw_x = handle["raw/X"]
        indptr = raw_x["indptr"][:]
        indices = raw_x["indices"]
        data = raw_x["data"]
        raw_valid_value_sum = 0.0
        total_blocks = (n_genes + args.gene_block_size - 1) // args.gene_block_size

        for block_number, gene_start in enumerate(
            range(0, n_genes, args.gene_block_size), start=1
        ):
            gene_end = min(gene_start + args.gene_block_size, n_genes)
            data_start = int(indptr[gene_start])
            data_end = int(indptr[gene_end])
            block_rows = indices[data_start:data_end]
            block_values = data[data_start:data_end]
            block_group = cell_to_group[block_rows]
            per_gene_lengths = np.diff(indptr[gene_start : gene_end + 1]).astype(
                np.int64
            )
            local_gene = np.repeat(
                np.arange(gene_end - gene_start, dtype=np.int64),
                per_gene_lengths,
            )
            valid = block_group >= 0
            raw_valid_value_sum += float(block_values[valid].sum())
            combined = (
                local_gene[valid] * n_groups + block_group[valid].astype(np.int64)
            )
            aggregated = np.bincount(
                combined,
                weights=block_values[valid],
                minlength=(gene_end - gene_start) * n_groups,
            ).reshape(gene_end - gene_start, n_groups)
            if not np.allclose(aggregated, np.round(aggregated)):
                raise ValueError(f"Non-integer aggregate detected at genes {gene_start}:{gene_end}")
            if aggregated.max(initial=0) > np.iinfo(np.int32).max:
                raise OverflowError("Pseudobulk count exceeds int32 range")
            all_counts[gene_start:gene_end, :] = np.rint(aggregated).astype(np.int32)
            if block_number == 1 or block_number % 25 == 0 or block_number == total_blocks:
                print(
                    f"AGGREGATION_PROGRESS block={block_number}/{total_blocks} "
                    f"genes={gene_end}/{n_genes}",
                    flush=True,
                )
        all_counts.flush()

    metadata["library_size"] = np.asarray(all_counts.sum(axis=0), dtype=np.int64)
    metadata.to_csv(out_dir / "all_group_metadata.csv", index=False)

    manifest_rows = []
    for cell_type in sorted(metadata["cell_type"].unique()):
        group_indices = metadata.index[metadata["cell_type"] == cell_type].to_numpy()
        cell_metadata = metadata.loc[group_indices].copy().reset_index(drop=True)
        slug = slugify(cell_type)
        binary_path = count_dir / f"{slug}.gene_by_sample.int32.bin"
        cell_counts = np.asarray(all_counts[:, group_indices], dtype=np.int32)
        cell_counts.T.tofile(binary_path)
        cell_metadata.to_csv(metadata_dir / f"{slug}.csv", index=False)
        manifest_rows.append(
            {
                "cell_type": cell_type,
                "slug": slug,
                "genes": len(gene_symbols),
                "samples": len(cell_metadata),
                "donors": cell_metadata["donor_id"].nunique(),
                "et_donors": cell_metadata.loc[
                    cell_metadata["disease"] == "essential tremor", "donor_id"
                ].nunique(),
                "normal_donors": cell_metadata.loc[
                    cell_metadata["disease"] == "normal", "donor_id"
                ].nunique(),
                "total_cells": int(cell_metadata["cell_count"].sum()),
                "total_counts": int(cell_counts.sum(dtype=np.int64)),
                "binary_bytes": binary_path.stat().st_size,
            }
        )
    manifest = pd.DataFrame(manifest_rows)
    manifest.to_csv(out_dir / "aggregation_manifest.csv", index=False)

    aggregated_count_sum = int(manifest["total_counts"].sum())
    qc = {
        "h5ad": str(h5ad_path),
        "minimum_cells_per_pseudobulk": args.min_cells,
        "genes": int(len(gene_symbols)),
        "pseudobulk_groups": int(len(metadata)),
        "cell_types": int(metadata["cell_type"].nunique()),
        "raw_valid_value_sum": int(round(raw_valid_value_sum)),
        "aggregated_count_sum": aggregated_count_sum,
        "count_conservation_pass": bool(
            int(round(raw_valid_value_sum)) == aggregated_count_sum
        ),
        "all_groups_binary": str(all_counts_path),
        "all_groups_binary_bytes": all_counts_path.stat().st_size,
        "binary_layout": "int32 little-endian; all-groups file is C-order gene-by-group; per-cell-type files are serialized for R column-major gene-by-sample matrices",
    }
    with open(out_dir / "aggregation_qc.json", "w", encoding="utf-8") as handle:
        json.dump(qc, handle, ensure_ascii=False, indent=2)
    print(json.dumps(qc, ensure_ascii=False), flush=True)


if __name__ == "__main__":
    main()
