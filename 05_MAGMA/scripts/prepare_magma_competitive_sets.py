from pathlib import Path

import pandas as pd


PROJECT_ROOT = Path(r"D:\CodexProjects\ET_MR_Stage0_20260821_v1")
VERSION_ROOT = PROJECT_ROOT / "outputs" / "v11_document_guided_strengthening_20260826"
RESULT_ROOT = VERSION_ROOT / "02_competitive_gene_set"

membership = pd.read_csv(RESULT_ROOT / "magma_gene_set_membership.csv", encoding="utf-8-sig")
genes = pd.read_csv(
    PROJECT_ROOT / "outputs" / "v4_final_evidence_20260822" / "22_magma_all_genes_35up10down.csv",
    encoding="utf-8-sig",
)[["gene_symbol", "magma_entrez_id"]].drop_duplicates("gene_symbol")

mapped = membership.merge(genes, on="gene_symbol", how="left")
mapped = mapped.dropna(subset=["magma_entrez_id"]).copy()
mapped["magma_entrez_id"] = mapped["magma_entrez_id"].astype(int)

set_path = RESULT_ROOT / "magma_competitive_sets.txt"
with set_path.open("w", encoding="utf-8", newline="\n") as handle:
    for set_name, group in mapped.groupby("gene_set", sort=False):
        identifiers = " ".join(str(value) for value in dict.fromkeys(group["magma_entrez_id"]))
        handle.write(f"{set_name} {identifiers}\n")

mapped.to_csv(
    RESULT_ROOT / "magma_competitive_set_membership.csv",
    index=False,
    encoding="utf-8-sig",
)
print(set_path)
