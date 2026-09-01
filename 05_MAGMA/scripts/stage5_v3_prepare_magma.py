from __future__ import annotations

import csv
import math
from pathlib import Path


ROOT = Path(r"D:\CodexProjects\ET_MR_Stage0_20260821_v1")
GWAS = ROOT / "raw" / "et_gwas" / "extracted_v1" / "G250_Essential_tremor_summary"
TARGETS = ROOT / "outputs" / "stage3e_prediction_target_union.csv"
GENE_LOC = ROOT / "tools" / "magma_v1.10" / "NCBI37.3.gene.loc"
WORK = ROOT / "work" / "v3_parallel_evidence_20260822" / "magma"
OUT = ROOT / "outputs" / "v3_parallel_evidence_20260822"
WORK.mkdir(parents=True, exist_ok=True)
OUT.mkdir(parents=True, exist_ok=True)


def prepare_pvalues() -> dict[str, int]:
    out_path = WORK / "et_2024_gwas_magma_pval.tsv"
    counts = {
        "input_rows": 0,
        "written_autosomal_rsid_rows": 0,
        "skipped_non_rsid": 0,
        "skipped_non_autosomal": 0,
        "skipped_invalid_p": 0,
    }
    with GWAS.open("r", encoding="utf-8", newline="") as src, out_path.open(
        "w", encoding="utf-8", newline=""
    ) as dst:
        header = src.readline().rstrip("\r\n").split()
        index = {name: i for i, name in enumerate(header)}
        dst.write("SNP\tP\n")
        for line in src:
            counts["input_rows"] += 1
            fields = line.split()
            rsid = fields[index["rsID"]]
            if not rsid.startswith("rs"):
                counts["skipped_non_rsid"] += 1
                continue
            chrom = fields[index["CHR"]].removeprefix("chr")
            if not chrom.isdigit() or not 1 <= int(chrom) <= 22:
                counts["skipped_non_autosomal"] += 1
                continue
            try:
                p_value = float(fields[index["P"]])
            except ValueError:
                counts["skipped_invalid_p"] += 1
                continue
            if not math.isfinite(p_value) or not 0 <= p_value <= 1:
                counts["skipped_invalid_p"] += 1
                continue
            dst.write(f"{rsid}\t{p_value:.16g}\n")
            counts["written_autosomal_rsid_rows"] += 1
    return counts


def prepare_target_set() -> dict[str, int]:
    symbol_to_ids: dict[str, list[str]] = {}
    gene_loc_lines: dict[str, str] = {}
    with GENE_LOC.open("r", encoding="utf-8", newline="") as handle:
        for line in handle:
            fields = line.rstrip("\r\n").split("\t")
            if len(fields) < 6:
                continue
            entrez, symbol = fields[0], fields[5]
            symbol_to_ids.setdefault(symbol.upper(), []).append(entrez)
            gene_loc_lines[entrez] = line.rstrip("\r\n")

    with TARGETS.open("r", encoding="utf-8-sig", newline="") as handle:
        targets = list(csv.DictReader(handle))

    mapped_ids: list[str] = []
    map_rows: list[dict[str, str]] = []
    for row in targets:
        symbol = row["gene_symbol"].strip()
        ids = symbol_to_ids.get(symbol.upper(), [])
        map_rows.append(
            {
                "gene_symbol": symbol,
                "magma_entrez_ids": ";".join(ids),
                "magma_symbol_mapping_status": "mapped" if ids else "not_mapped",
            }
        )
        mapped_ids.extend(ids)

    mapped_ids = list(dict.fromkeys(mapped_ids))
    with (WORK / "predicted_105_magma_set.txt").open("w", encoding="utf-8") as handle:
        handle.write("ET_PREDICTED_TARGET_UNION_105 " + " ".join(mapped_ids) + "\n")
    with (WORK / "predicted_105_NCBI37.3.gene.loc").open("w", encoding="utf-8") as handle:
        for entrez in mapped_ids:
            handle.write(gene_loc_lines[entrez] + "\n")

    with (OUT / "11_predicted_105_magma_id_mapping.csv").open(
        "w", encoding="utf-8-sig", newline=""
    ) as handle:
        writer = csv.DictWriter(handle, fieldnames=list(map_rows[0]))
        writer.writeheader()
        writer.writerows(map_rows)

    return {
        "target_symbols": len(targets),
        "mapped_target_symbols": sum(
            row["magma_symbol_mapping_status"] == "mapped" for row in map_rows
        ),
        "unique_entrez_ids": len(mapped_ids),
    }


def main() -> None:
    p_counts = prepare_pvalues()
    target_counts = prepare_target_set()
    qa_path = OUT / "12_magma_input_qa.csv"
    rows = [
        {"metric": key, "value": value}
        for key, value in {**p_counts, **target_counts}.items()
    ]
    with qa_path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["metric", "value"])
        writer.writeheader()
        writer.writerows(rows)
    print({**p_counts, **target_counts})


if __name__ == "__main__":
    main()
