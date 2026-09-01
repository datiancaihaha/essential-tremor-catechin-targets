from pathlib import Path

import numpy as np
import pandas as pd


PROJECT_ROOT = Path(r"D:\CodexProjects\ET_MR_Stage0_20260821_v1")
VERSION_ROOT = PROJECT_ROOT / "outputs" / "v11_document_guided_strengthening_20260826"
RESULT_ROOT = VERSION_ROOT / "02_competitive_gene_set"


def bh_adjust(values: pd.Series) -> pd.Series:
    values = values.astype(float)
    order = np.argsort(values.to_numpy())
    ranked = values.to_numpy()[order]
    adjusted = np.minimum.accumulate(
        (ranked * len(ranked) / np.arange(1, len(ranked) + 1))[::-1]
    )[::-1]
    adjusted = np.minimum(adjusted, 1.0)
    inverse = np.empty_like(order)
    inverse[order] = np.arange(len(order))
    return pd.Series(adjusted[inverse], index=values.index)


result_path = RESULT_ROOT / "magma_competitive.gsa.out.txt"
rows = []
with result_path.open("r", encoding="utf-8") as handle:
    for line in handle:
        if not line.strip() or line.startswith("#") or line.startswith("VARIABLE"):
            continue
        fields = line.split()
        rows.append(
            {
                "gene_set": fields[-1],
                "test_type": fields[1],
                "tested_gene_count": int(fields[2]),
                "beta": float(fields[3]),
                "standardized_beta": float(fields[4]),
                "standard_error": float(fields[5]),
                "pvalue": float(fields[6]),
            }
        )

results = pd.DataFrame(rows)
results["nominal_p_lt_0_05"] = results["pvalue"] < 0.05
results["bh_fdr_supplementary"] = bh_adjust(results["pvalue"])
results["candidate_retention_rule"] = "nominal P < 0.05"
results["multiplicity_role"] = "supplementary evidence strength only"
results.to_csv(
    RESULT_ROOT / "magma_competitive_gene_set_results.csv",
    index=False,
    encoding="utf-8-sig",
)
print(results.to_string(index=False))
