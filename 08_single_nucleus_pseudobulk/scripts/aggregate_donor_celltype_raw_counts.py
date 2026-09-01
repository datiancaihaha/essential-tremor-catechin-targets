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
    parser.add_argument("--cell-counts", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--gene-block-size", type=int, default=64)
    args = parser.parse_args()

    h5ad_path = Path(args.h5ad)
    out_dir = Path(args.out)
    count_dir = out_dir / "count_matrices"
    metadata_dir = out_dir / "sample_metadata"
    for directory in (out_dir, count_dir, metadata_dir):
        directory.mkdir(parents=True, exist_ok=True)

    batch_counts = pd.read_csv(
        args.cell_counts,
        dtype={"donor_id": str, "sample_id": str},
    )
    required = {
        "donor_id",
        "cell_type",
        "seq_batch",
        "cell_count",
        "disease",
        "sex",
        "age",
    }
    if not required.issubset(batch_counts.columns):
        raise ValueError("The cell-count table is missing required columns")
    if (batch_counts["cell_count"] < 0).any():
        raise ValueError("Cell counts must be non-negative")

    donor_metadata = batch_counts[
        ["donor_id", "disease", "sex", "age"]
    ].drop_duplicates()
    if donor_metadata["donor_id"].duplicated().any():
        raise ValueError("Donor metadata are not unique")
    batch_profiles = (
        batch_counts[["donor_id", "seq_batch"]]
        .drop_duplicates()
        .groupby("donor_id", sort=False)["seq_batch"]
        .agg(lambda values: "+".join(sorted(set(values))))
        .rename("batch_profile")
        .reset_index()
    )
    donor_metadata = donor_metadata.merge(
        batch_profiles,
        on="donor_id",
        how="left",
        validate="one_to_one",
    )
    donor_cell_counts = (
        batch_counts.groupby(
            ["donor_id", "cell_type"], as_index=False, sort=True
        )["cell_count"]
        .sum()
        .merge(donor_metadata, on="donor_id", how="left", validate="many_to_one")
    )
    donor_cell_counts = donor_cell_counts.sort_values(
        ["cell_type", "donor_id"], kind="stable"
    ).reset_index(drop=True)
    donor_cell_counts["group_index"] = np.arange(
        len(donor_cell_counts), dtype=np.int32
    )

    with h5py.File(h5ad_path, "r") as handle:
        donor_names, donor_codes = categorical(handle["obs/donor_id"])
        cell_type_names, cell_type_codes = categorical(handle["obs/cell_type"])
        donor_lookup = {value: index for index, value in enumerate(donor_names)}
        cell_type_lookup = {
            value: index for index, value in enumerate(cell_type_names)
        }
        if set(donor_cell_counts["donor_id"]) != set(donor_names):
            raise ValueError("Donor identifiers differ between metadata and H5AD")
        if set(donor_cell_counts["cell_type"]) != set(cell_type_names):
            raise ValueError("Cell-type identifiers differ between metadata and H5AD")

        n_full_groups = len(donor_names) * len(cell_type_names)
        full_to_compact = np.full(n_full_groups, -1, dtype=np.int32)
        for row in donor_cell_counts.itertuples(index=False):
            full_code = (
                donor_lookup[row.donor_id] * len(cell_type_names)
                + cell_type_lookup[row.cell_type]
            )
            full_to_compact[full_code] = row.group_index
        cell_full_groups = donor_codes * len(cell_type_names) + cell_type_codes
        cell_to_group = full_to_compact[cell_full_groups]
        if (cell_to_group < 0).any():
            raise ValueError("At least one nucleus could not be assigned to a donor-cell-type group")

        observed_cells = np.bincount(
            cell_to_group,
            minlength=len(donor_cell_counts),
        )
        if not np.array_equal(
            observed_cells.astype(np.int64),
            donor_cell_counts["cell_count"].to_numpy(dtype=np.int64),
        ):
            raise ValueError("Donor-cell-type nucleus counts do not match the H5AD annotations")

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
        n_groups = len(donor_cell_counts)
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
        raw_value_sum = 0.0
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
            per_gene_lengths = np.diff(
                indptr[gene_start : gene_end + 1]
            ).astype(np.int64)
            local_gene = np.repeat(
                np.arange(gene_end - gene_start, dtype=np.int64),
                per_gene_lengths,
            )
            raw_value_sum += float(block_values.sum())
            combined = (
                local_gene * n_groups + block_group.astype(np.int64)
            )
            aggregated = np.bincount(
                combined,
                weights=block_values,
                minlength=(gene_end - gene_start) * n_groups,
            ).reshape(gene_end - gene_start, n_groups)
            if not np.allclose(aggregated, np.round(aggregated)):
                raise ValueError(
                    f"Non-integer aggregate at genes {gene_start}:{gene_end}"
                )
            if aggregated.max(initial=0) > np.iinfo(np.int32).max:
                raise OverflowError("A donor-cell-type count exceeds int32 range")
            all_counts[gene_start:gene_end, :] = np.rint(aggregated).astype(np.int32)
            if (
                block_number == 1
                or block_number % 25 == 0
                or block_number == total_blocks
            ):
                print(
                    f"AGGREGATION_PROGRESS block={block_number}/{total_blocks} "
                    f"genes={gene_end}/{n_genes}",
                    flush=True,
                )
        all_counts.flush()

    donor_cell_counts["library_size"] = np.asarray(
        all_counts.sum(axis=0), dtype=np.int64
    )
    donor_cell_counts.to_csv(out_dir / "all_group_metadata.csv", index=False)

    manifest_rows = []
    for cell_type in sorted(donor_cell_counts["cell_type"].unique()):
        group_indices = donor_cell_counts.index[
            donor_cell_counts["cell_type"] == cell_type
        ].to_numpy()
        cell_metadata = donor_cell_counts.loc[group_indices].copy().reset_index(drop=True)
        slug = slugify(cell_type)
        binary_path = count_dir / f"{slug}.gene_by_donor.int32.bin"
        cell_counts = np.asarray(all_counts[:, group_indices], dtype=np.int32)
        cell_counts.T.tofile(binary_path)
        cell_metadata.to_csv(metadata_dir / f"{slug}.csv", index=False)
        manifest_rows.append(
            {
                "cell_type": cell_type,
                "slug": slug,
                "genes": len(gene_symbols),
                "donors": len(cell_metadata),
                "et_donors": int(
                    (cell_metadata["disease"] == "essential tremor").sum()
                ),
                "control_donors": int((cell_metadata["disease"] == "normal").sum()),
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
        "aggregation_unit": "donor by author-defined cell type",
        "genes": int(len(gene_symbols)),
        "donor_cell_type_groups": int(len(donor_cell_counts)),
        "donors": int(donor_cell_counts["donor_id"].nunique()),
        "cell_types": int(donor_cell_counts["cell_type"].nunique()),
        "raw_value_sum": int(round(raw_value_sum)),
        "aggregated_count_sum": aggregated_count_sum,
        "count_conservation_pass": bool(
            int(round(raw_value_sum)) == aggregated_count_sum
        ),
        "cell_count_annotation_match_pass": True,
        "all_groups_binary": str(all_counts_path),
        "all_groups_binary_bytes": all_counts_path.stat().st_size,
        "binary_layout": (
            "int32 little-endian; all-groups file is C-order gene-by-group; "
            "per-cell-type files are serialized for R column-major gene-by-donor matrices"
        ),
    }
    with open(out_dir / "aggregation_qc.json", "w", encoding="utf-8") as handle:
        json.dump(qc, handle, ensure_ascii=False, indent=2)
    print(json.dumps(qc, ensure_ascii=False), flush=True)


if __name__ == "__main__":
    main()
