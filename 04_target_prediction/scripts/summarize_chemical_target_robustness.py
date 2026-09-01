from __future__ import annotations

import csv
from collections import defaultdict
from pathlib import Path


PROJECT_ROOT = Path(r"D:\CodexProjects\ET_MR_Stage0_20260821_v1")
RESULT_ROOT = (
    PROJECT_ROOT
    / "outputs"
    / "v11_document_guided_strengthening_20260826"
    / "03_chemical_forms"
)


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def valid_gene_symbol(value: str) -> bool:
    return value.strip().upper() not in {"", "N/A", "NA", "-"}


classification = read_csv(RESULT_ROOT / "chemical_exposure_classification.csv")
predictions = read_csv(RESULT_ROOT / "swisstargetprediction_all_chemical_forms.csv")
experimental = read_csv(RESULT_ROOT / "experimental_target_database_summary.csv")
activities = read_csv(RESULT_ROOT / "chembl_experimental_activities.csv")

form_order = [row["record_id"] for row in classification]
grade_a_forms = [
    row["record_id"] for row in classification if row["family_grade"].strip() == "A"
]
if len(grade_a_forms) != 4:
    raise RuntimeError(f"Expected four Grade-A circulating conjugate forms; found {grade_a_forms}")

prediction_by_form_gene: dict[tuple[str, str], dict[str, str]] = {}
genes: set[str] = set()
for row in predictions:
    gene = row["gene_symbol"].strip()
    if not valid_gene_symbol(gene):
        continue
    genes.add(gene)
    prediction_by_form_gene[(row["record_id"], gene)] = row

matrix_rows: list[dict[str, object]] = []
for gene in sorted(genes):
    row: dict[str, object] = {"gene_symbol": gene}
    for form in form_order:
        hit = prediction_by_form_gene.get((form, gene))
        row[f"{form}_predicted"] = bool(hit)
        row[f"{form}_rank"] = int(hit["rank"]) if hit else ""
        row[f"{form}_probability"] = float(hit["probability"]) if hit else ""
    matrix_rows.append(row)

matrix_fields = ["gene_symbol"]
for form in form_order:
    matrix_fields.extend(
        [f"{form}_predicted", f"{form}_rank", f"{form}_probability"]
    )
write_csv(
    RESULT_ROOT / "exact_form_target_prediction_matrix.csv",
    matrix_rows,
    matrix_fields,
)

all_form_gene_sets = {
    form: {
        row["gene_symbol"].strip()
        for row in predictions
        if row["record_id"] == form and valid_gene_symbol(row["gene_symbol"])
    }
    for form in form_order
}
all_form_union = set().union(*all_form_gene_sets.values())
all_form_intersection = set.intersection(*all_form_gene_sets.values())
write_csv(
    RESULT_ROOT / "chemical_form_target_summary.csv",
    [
        {
            "record_id": form,
            "valid_unique_gene_count": len(all_form_gene_sets[form]),
            "CA1_predicted": "CA1" in all_form_gene_sets[form],
            "CA2_predicted": "CA2" in all_form_gene_sets[form],
            "CA3_predicted": "CA3" in all_form_gene_sets[form],
            "CA9_predicted": "CA9" in all_form_gene_sets[form],
            "CA12_predicted": "CA12" in all_form_gene_sets[form],
        }
        for form in form_order
    ],
    [
        "record_id",
        "valid_unique_gene_count",
        "CA1_predicted",
        "CA2_predicted",
        "CA3_predicted",
        "CA9_predicted",
        "CA12_predicted",
    ],
)
write_csv(
    RESULT_ROOT / "all_chemical_form_target_overlap_summary.csv",
    [
        {
            "form_count": len(form_order),
            "target_union_count": len(all_form_union),
            "target_intersection_count": len(all_form_intersection),
            "intersection_genes": ";".join(sorted(all_form_intersection)),
        }
    ],
    [
        "form_count",
        "target_union_count",
        "target_intersection_count",
        "intersection_genes",
    ],
)

grade_a_rows: list[dict[str, object]] = []
grade_a_genes = sorted(
    {
        row["gene_symbol"].strip()
        for row in predictions
        if row["record_id"] in grade_a_forms
        and valid_gene_symbol(row["gene_symbol"])
    }
)
for gene in grade_a_genes:
    hits = [prediction_by_form_gene.get((form, gene)) for form in grade_a_forms]
    observed = [hit for hit in hits if hit is not None]
    ranks = [int(hit["rank"]) for hit in observed]
    probabilities = [float(hit["probability"]) for hit in observed]
    present = {form: prediction_by_form_gene.get((form, gene)) is not None for form in grade_a_forms}
    grade_a_rows.append(
        {
            "gene_symbol": gene,
            "predicted_form_count": len(observed),
            "predicted_in_all_four_forms": len(observed) == 4,
            "predicted_in_both_sulfates": present["LACTONE_3_SULFATE"]
            and present["LACTONE_4_SULFATE"],
            "predicted_in_both_glucuronides": present["LACTONE_3_GLUCURONIDE"]
            and present["LACTONE_4_GLUCURONIDE"],
            "best_rank": min(ranks),
            "worst_rank": max(ranks),
            "maximum_probability": max(probabilities),
            "mean_probability_across_present_forms": sum(probabilities) / len(probabilities),
        }
    )

grade_a_rows.sort(
    key=lambda row: (
        -int(row["predicted_form_count"]),
        int(row["best_rank"]),
        str(row["gene_symbol"]),
    )
)
write_csv(
    RESULT_ROOT / "circulating_conjugate_target_robustness.csv",
    grade_a_rows,
    [
        "gene_symbol",
        "predicted_form_count",
        "predicted_in_all_four_forms",
        "predicted_in_both_sulfates",
        "predicted_in_both_glucuronides",
        "best_rank",
        "worst_rank",
        "maximum_probability",
        "mean_probability_across_present_forms",
    ],
)

family_groups = {
    "All circulating conjugates": grade_a_forms,
    "Sulfate conjugates": ["LACTONE_3_SULFATE", "LACTONE_4_SULFATE"],
    "Glucuronide conjugates": ["LACTONE_3_GLUCURONIDE", "LACTONE_4_GLUCURONIDE"],
}
family_rows: list[dict[str, object]] = []
for family, forms in family_groups.items():
    sets = [
        {
            row["gene_symbol"].strip()
            for row in predictions
            if row["record_id"] == form and valid_gene_symbol(row["gene_symbol"])
        }
        for form in forms
    ]
    union = set().union(*sets)
    intersection = set.intersection(*sets)
    family_rows.append(
        {
            "chemical_family": family,
            "form_count": len(forms),
            "target_union_count": len(union),
            "target_intersection_count": len(intersection),
            "intersection_genes": ";".join(sorted(intersection)),
        }
    )
write_csv(
    RESULT_ROOT / "chemical_family_target_overlap_summary.csv",
    family_rows,
    [
        "chemical_family",
        "form_count",
        "target_union_count",
        "target_intersection_count",
        "intersection_genes",
    ],
)

activity_by_record: dict[str, list[dict[str, str]]] = defaultdict(list)
for row in activities:
    activity_by_record[row["record_id"]].append(row)

evidence_rows: list[dict[str, object]] = []
for row in experimental:
    record_activities = activity_by_record[row["record_id"]]
    target_resolved = [
        activity
        for activity in record_activities
        if activity["target_pref_name"].strip()
        and activity["target_pref_name"].strip().lower() not in {"unchecked", "ht-29"}
    ]
    evidence_rows.append(
        {
            **row,
            "chembl_target_resolved_activity_count": len(target_resolved),
            "target_resolved_experimental_evidence": len(target_resolved) > 0,
            "evidence_interpretation": (
                "Target-resolved experimental activity available"
                if target_resolved
                else "No target-resolved experimental activity identified"
            ),
        }
    )
write_csv(
    RESULT_ROOT / "experimental_target_evidence_classification.csv",
    evidence_rows,
    list(experimental[0].keys())
    + [
        "chembl_target_resolved_activity_count",
        "target_resolved_experimental_evidence",
        "evidence_interpretation",
    ],
)

print(
    "Chemical-form target robustness completed: "
    f"{len(grade_a_genes)} Grade-A family-union targets; "
    f"{sum(row['predicted_in_all_four_forms'] for row in grade_a_rows)} present in all four forms."
)
