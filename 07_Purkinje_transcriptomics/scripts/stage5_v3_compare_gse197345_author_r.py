from __future__ import annotations

import csv
import math
import xml.etree.ElementTree as ET
from pathlib import Path


ROOT = Path(r"D:\CodexProjects\ET_MR_Stage0_20260821_v1")
BIOC = (
    ROOT
    / "raw"
    / "stage5_cerebellum_validation_v2_20260822"
    / "literature"
    / "PMC10359949.bioc.xml"
)
R_RESULTS = ROOT / "outputs" / "v3_parallel_evidence_20260822" / "08_gse197345_r_deseq2_all.csv"
TARGETS = ROOT / "outputs" / "stage3e_prediction_target_union.csv"
OUT = ROOT / "outputs" / "v3_parallel_evidence_20260822"
OUT.mkdir(parents=True, exist_ok=True)


def clean_number(value: str) -> float:
    return float(value.replace("−", "-").replace("–", "-").strip())


def read_author_table() -> list[dict[str, object]]:
    document = ET.parse(BIOC).getroot()
    table_xml = None
    for passage in document.iter("passage"):
        info = {node.attrib.get("key"): "".join(node.itertext()) for node in passage.findall("infon")}
        if info.get("id") == "T1" and info.get("type") == "table":
            table_xml = info.get("xml")
            break
    if not table_xml:
        raise RuntimeError("Could not find author Table 1 in PMC10359949 BioC XML")

    table = ET.fromstring(table_xml)
    rows = []
    for tr in table.findall("./tbody/tr"):
        cells = ["".join(td.itertext()).strip() for td in tr.findall("td")]
        if len(cells) != 6:
            raise RuntimeError(f"Unexpected Table 1 row: {cells}")
        rows.append(
            {
                "gene_symbol": cells[0],
                "gene_name": cells[1],
                "author_baseMean": clean_number(cells[2]),
                "author_log2FoldChange": clean_number(cells[3]),
                "author_pvalue": clean_number(cells[4]),
                "author_padj": clean_number(cells[5]),
                "author_result_scope": "published_Table_1_BH_FDR_lt_0_05_only",
            }
        )
    if len(rows) != 36:
        raise RuntimeError(f"Expected 36 author rows, found {len(rows)}")
    return rows


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def as_float(value: str) -> float:
    if value == "" or value.lower() in {"nan", "na"}:
        return math.nan
    return float(value)


def main() -> None:
    author = read_author_table()
    write_csv(OUT / "13_gse197345_author_table1_36.csv", author)

    r_rows = read_csv(R_RESULTS)
    r_by_gene = {row["gene_symbol"]: row for row in r_rows}
    comparison = []
    for row in author:
        r = r_by_gene.get(str(row["gene_symbol"]))
        comparison.append(
            {
                **row,
                "r_baseMean": as_float(r["baseMean"]) if r else math.nan,
                "r_log2FoldChange": as_float(r["log2FoldChange"]) if r else math.nan,
                "r_pvalue": as_float(r["pvalue"]) if r else math.nan,
                "r_padj": as_float(r["padj"]) if r else math.nan,
                "r_bh_fdr_lt_0_05": bool(r and r["bh_fdr_lt_0_05_supplementary"] == "TRUE"),
            }
        )
    write_csv(OUT / "14_gse197345_author_vs_r_deseq2_36.csv", comparison)

    targets = read_csv(TARGETS)
    target_rows = []
    author_genes = {str(row["gene_symbol"]) for row in author}
    for target in targets:
        gene = target["gene_symbol"]
        r = r_by_gene.get(gene)
        target_rows.append(
            {
                "gene_symbol": gene,
                "in_swiss": target["in_swiss"],
                "in_sea": target["in_sea"],
                "prediction_sources": target["prediction_sources"],
                "purkinje_r_baseMean": as_float(r["baseMean"]) if r else math.nan,
                "purkinje_r_log2FoldChange": as_float(r["log2FoldChange"]) if r else math.nan,
                "purkinje_r_pvalue": as_float(r["pvalue"]) if r else math.nan,
                "purkinje_r_padj_supplementary": as_float(r["padj"]) if r else math.nan,
                "purkinje_r_nominal_p_lt_0_05": bool(r and r["nominal_p_lt_0_05"] == "TRUE"),
                "in_author_fdr_table1": gene in author_genes,
                "author_full_result_availability": "FDR-significant rows only; complete author result table not public",
            }
        )
    write_csv(OUT / "15_gse197345_prediction105_r_deseq2.csv", target_rows)

    r_fdr_genes = {
        row["gene_symbol"] for row in r_rows if row["bh_fdr_lt_0_05_supplementary"] == "TRUE"
    }
    author_genes = {str(row["gene_symbol"]) for row in author}
    qa = [
        {"metric": "author_table1_fdr_genes", "value": len(author_genes)},
        {"metric": "r_deseq2_fdr_genes", "value": len(r_fdr_genes)},
        {"metric": "author_genes_recovered_by_r_fdr", "value": len(author_genes & r_fdr_genes)},
        {"metric": "r_only_fdr_genes", "value": len(r_fdr_genes - author_genes)},
        {"metric": "author_only_fdr_genes", "value": len(author_genes - r_fdr_genes)},
        {
            "metric": "prediction105_r_nominal_p_lt_0_05",
            "value": sum(bool(row["purkinje_r_nominal_p_lt_0_05"]) for row in target_rows),
        },
        {
            "metric": "prediction105_in_author_fdr_table1",
            "value": sum(bool(row["in_author_fdr_table1"]) for row in target_rows),
        },
    ]
    write_csv(OUT / "16_gse197345_author_r_comparison_qa.csv", qa)
    print(qa)


if __name__ == "__main__":
    main()
