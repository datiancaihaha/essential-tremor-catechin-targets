from __future__ import annotations

import importlib.util
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd


def load_stage0_helpers(root: Path):
    path = root / "scripts" / "stage0_mibiogen_gutmgene_audit.py"
    spec = importlib.util.spec_from_file_location("stage0_mibiogen_gutmgene_audit", path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Cannot load stage0 helpers from {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def main() -> None:
    root = Path(sys.argv[1]).resolve() if len(sys.argv) > 1 else Path.cwd().resolve()
    output_dir = root / "outputs"
    qa_dir = root / "qa"
    qa_dir.mkdir(parents=True, exist_ok=True)

    ivw = pd.read_csv(output_dir / "stage2_ivw_primary_results.csv")
    methods = pd.read_csv(output_dir / "stage2_mr_results_all_methods.csv")
    sensitivity = pd.read_csv(output_dir / "stage2_sensitivity_summary.csv")
    loo_summary = pd.read_csv(output_dir / "stage2_leave_one_out_summary.csv")
    mr_presso = pd.read_csv(output_dir / "stage2_mr_presso_summary.csv")
    stage1_ivs = pd.read_csv(output_dir / "stage1_clumped_independent_ivs.csv")
    coverage = pd.read_csv(output_dir / "gutmgene_taxa_coverage.csv")
    strict_edges = pd.read_csv(output_dir / "gutmgene_strict_human_microbe_metabolite_edges.csv")

    selected = ivw.loc[ivw["p_value"] < 0.05].copy()
    if len(selected) != 10:
        raise RuntimeError(f"Expected 10 nominal IVW hits, found {len(selected)}")
    selected["selected_for_downstream"] = "yes"
    selected["selection_rule"] = "IVW nominal P<0.05; BH-FDR reported but not used as the selection gate"
    selected["selection_rule_zh"] = "按用户指定规则：IVW名义P<0.05入选；BH-FDR仅报告，不作为入选门槛"

    method_subset = methods.loc[
        methods["method"].isin(["Weighted median", "MR Egger", "Weighted mode"]),
        ["bac", "method", "beta", "se", "p_value", "odds_ratio", "or_ci_lower", "or_ci_upper", "method_status"],
    ].copy()
    method_subset["method_key"] = method_subset["method"].str.lower().str.replace(r"[^a-z0-9]+", "_", regex=True)
    wide_parts = []
    for value in ["beta", "se", "p_value", "odds_ratio", "or_ci_lower", "or_ci_upper", "method_status"]:
        pivot = method_subset.pivot(index="bac", columns="method_key", values=value)
        pivot.columns = [f"{method}_{value}" for method in pivot.columns]
        wide_parts.append(pivot)
    method_wide = pd.concat(wide_parts, axis=1).reset_index()

    selected = selected.merge(method_wide, on="bac", how="left", validate="one_to_one")
    selected = selected.merge(sensitivity, on="bac", how="left", validate="one_to_one", suffixes=("", "_sensitivity"))
    selected = selected.merge(loo_summary, on="bac", how="left", validate="one_to_one")
    selected = selected.merge(coverage, on="bac", how="left", validate="one_to_one")
    presso_columns = [
        "bac", "run_status", "run_status_zh", "global_p_text", "global_p_numeric",
        "global_test_significant", "outlier_count", "outlier_rsids", "corrected_beta", "corrected_p",
    ]
    presso_subset = mr_presso[presso_columns].rename(
        columns={column: f"mr_presso_{column}" for column in presso_columns if column != "bac"}
    )
    selected = selected.merge(presso_subset, on="bac", how="left", validate="one_to_one")

    sensitivity_methods = ["weighted_median", "mr_egger", "weighted_mode"]
    selected["direction_concordant_all_sensitivity_methods"] = True
    selected["nominal_supporting_sensitivity_method_count"] = 0
    for method in sensitivity_methods:
        same_direction = np.sign(selected[f"{method}_beta"]) == np.sign(selected["beta"])
        selected["direction_concordant_all_sensitivity_methods"] &= same_direction
        selected["nominal_supporting_sensitivity_method_count"] += (
            same_direction & (selected[f"{method}_p_value"] < 0.05)
        ).astype(int)
    selected["no_ivw_heterogeneity_signal"] = selected["ivw_q_p"] >= 0.05
    selected["no_egger_intercept_signal"] = selected["egger_intercept_p"] >= 0.05
    selected["no_mr_presso_global_signal"] = ~selected["mr_presso_global_test_significant"].eq(True)
    selected["loo_direction_stable"] = selected["all_loo_direction_concordant"].eq(True)
    diagnostic_clean = (
        selected["no_ivw_heterogeneity_signal"]
        & selected["no_egger_intercept_signal"]
        & selected["no_mr_presso_global_signal"]
        & selected["loo_direction_stable"]
    )
    selected["robustness_status"] = np.select(
        [
            diagnostic_clean
            & selected["direction_concordant_all_sensitivity_methods"]
            & (selected["nominal_supporting_sensitivity_method_count"] >= 1),
            diagnostic_clean & selected["direction_concordant_all_sensitivity_methods"],
            diagnostic_clean,
        ],
        [
            "NOMINAL_WITH_MULTI_METHOD_SUPPORT",
            "NOMINAL_DIRECTIONALLY_CONCORDANT",
            "NOMINAL_IVW_ONLY_DIRECTION_DISCORDANCE",
        ],
        default="NOMINAL_WITH_SENSITIVITY_WARNING",
    )
    robustness_zh = {
        "NOMINAL_WITH_MULTI_METHOD_SUPPORT": "名义入选；至少一种稳健方法同向且P<0.05，敏感性诊断未见警示",
        "NOMINAL_DIRECTIONALLY_CONCORDANT": "名义入选；三种敏感性方法方向一致，但均未达P<0.05，诊断未见警示",
        "NOMINAL_IVW_ONLY_DIRECTION_DISCORDANCE": "名义入选；至少一种敏感性方法方向不一致，诊断未见其他警示",
        "NOMINAL_WITH_SENSITIVITY_WARNING": "名义入选；存在异质性、多效性或留一法警示",
    }
    selected["robustness_status_zh"] = selected["robustness_status"].map(robustness_zh)

    selected_iv_rows = stage1_ivs.loc[
        stage1_ivs["mr_usable"].eq("yes") & stage1_ivs["bac"].isin(selected["bac"]),
        ["bac", "analysis_rsid"],
    ]
    iv_sets = {
        taxon: tuple(sorted(group["analysis_rsid"].astype(str)))
        for taxon, group in selected_iv_rows.groupby("bac", sort=True)
    }
    unique_sets = sorted(set(iv_sets.values()))
    iv_set_id = {value: f"SELECTED_IVSET_{index:03d}" for index, value in enumerate(unique_sets, start=1)}
    members_by_set: dict[tuple[str, ...], list[str]] = {}
    for taxon, value in iv_sets.items():
        members_by_set.setdefault(value, []).append(taxon)
    selected["selected_iv_set_id"] = selected["bac"].map(lambda taxon: iv_set_id[iv_sets[taxon]])
    selected["identical_selected_iv_set_group_size"] = selected["bac"].map(lambda taxon: len(members_by_set[iv_sets[taxon]]))
    selected["identical_selected_iv_set_members"] = selected["bac"].map(
        lambda taxon: ";".join(sorted(members_by_set[iv_sets[taxon]]))
    )
    selected["non_independent_identical_iv_set"] = np.where(
        selected["identical_selected_iv_set_group_size"] > 1, "yes", "no"
    )
    selected = selected.sort_values(["p_value", "bac"], kind="stable")

    selected_bac = set(selected["bac"])
    strict_selected = strict_edges.loc[strict_edges["mibiogen_bac"].isin(selected_bac)].copy()
    helpers = load_stage0_helpers(root)
    strict_selected["is_gaba"] = strict_selected["Metabolite"].map(helpers.is_gaba)
    strict_selected["mapping_scope"] = "STRICT_SAME_TAXONOMIC_RANK"
    strict_selected["mapping_scope_zh"] = "严格映射：同一分类层级且规范化名称完全匹配"
    strict_selected.to_csv(
        output_dir / "stage2_nominal_selected_gutmgene_strict_edges.csv",
        index=False,
        encoding="utf-8-sig",
    )

    raw_edges, raw_encoding = helpers.read_csv(root / "raw" / "gutmgene" / "Gut Microbe-Microbial metabolite.csv")
    selected_lookup = selected.set_index("bac")
    permissive_rows: list[dict] = []
    for row in raw_edges:
        if helpers.norm(row.get("human/mouse", "")) != "human":
            continue
        if helpers.norm(row.get("Rank", "")) != "species":
            continue
        microbe_name = row.get("Gut Microbiota", "")
        genus_key = helpers.norm(microbe_name).split(" ", 1)[0]
        for bac, selected_row in selected_lookup.iterrows():
            if selected_row.get("rank") != "genus":
                continue
            if int(selected_row.get("strict_human_microbe_metabolite_rows", 0)) > 0:
                continue
            if genus_key != helpers.norm(selected_row.get("taxon_name", "")):
                continue
            permissive_rows.append(
                {
                    "selected_mibiogen_bac": bac,
                    "selected_genus": selected_row.get("taxon_name", ""),
                    "mapping_scope": "PERMISSIVE_GENUS_TO_SPECIES",
                    "mapping_scope_zh": "宽松映射：MiBioGen属名对应gutMGene物种首词，仅作探索性上界",
                    "is_gaba": helpers.is_gaba(row.get("Metabolite", "")),
                    **row,
                }
            )
    permissive = pd.DataFrame(permissive_rows)
    if not permissive.empty:
        permissive = permissive.sort_values(
            ["selected_mibiogen_bac", "Gut Microbiota", "Metabolite", "PMID"], kind="stable"
        )
    permissive.to_csv(
        output_dir / "stage2_nominal_selected_gutmgene_permissive_species_edges.csv",
        index=False,
        encoding="utf-8-sig",
    )

    strict_non_gaba = strict_selected.loc[~strict_selected["is_gaba"]]
    permissive_non_gaba = permissive.loc[~permissive["is_gaba"]] if not permissive.empty else permissive
    strict_non_gaba_bac = set(strict_non_gaba["mibiogen_bac"])
    permissive_non_gaba_bac = set(permissive_non_gaba["selected_mibiogen_bac"]) if not permissive_non_gaba.empty else set()
    selected["gutmgene_pathway_priority"] = np.select(
        [
            selected["bac"].isin(strict_non_gaba_bac),
            selected["bac"].isin(permissive_non_gaba_bac),
        ],
        ["STRICT_NON_GABA_PRIORITY", "PERMISSIVE_SPECIES_EXPLORATORY"],
        default="NO_MATCHED_GUTMGENE_METABOLITE_EDGE",
    )
    priority_zh = {
        "STRICT_NON_GABA_PRIORITY": "严格同层级匹配到非GABA代谢物，可进入优先通路核查",
        "PERMISSIVE_SPECIES_EXPLORATORY": "仅属到物种宽松匹配到非GABA代谢物，只作探索性候选",
        "NO_MATCHED_GUTMGENE_METABOLITE_EDGE": "未匹配到可用的人源gutMGene微生物-代谢物边",
    }
    selected["gutmgene_pathway_priority_zh"] = selected["gutmgene_pathway_priority"].map(priority_zh)
    selected.to_csv(output_dir / "stage2_nominal_selected_taxa.csv", index=False, encoding="utf-8-sig")
    metrics = {
        "status": "PASS_NOMINAL_SELECTION_WITH_LIMITED_STRICT_GUTMGENE_COVERAGE",
        "status_zh": "通过名义筛选；但严格gutMGene覆盖有限，宽松物种映射仅供探索",
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "selection_rule": "IVW nominal P<0.05; BH-FDR reported but not used as the selection gate",
        "selection_rule_zh": "按用户指定规则：IVW名义P<0.05入选；BH-FDR仅报告，不作为入选门槛",
        "counts": {
            "selected_taxa": int(len(selected)),
            "selected_taxa_bh_fdr_lt_0_05": int((selected["fdr_bh"] < 0.05).sum()),
            "strict_covered_selected_taxa": int(selected["strict_covered"].eq("yes").sum()),
            "strict_selected_edges": int(len(strict_selected)),
            "strict_selected_non_gaba_edges": int(len(strict_non_gaba)),
            "strict_selected_unique_non_gaba_metabolites": int(strict_non_gaba["Metabolite"].nunique()),
            "permissive_species_edges": int(len(permissive)),
            "permissive_species_non_gaba_edges": int(len(permissive_non_gaba)),
            "permissive_selected_taxa": int(permissive["selected_mibiogen_bac"].nunique()) if not permissive.empty else 0,
            "nominal_with_multi_method_support": int(selected["robustness_status"].eq("NOMINAL_WITH_MULTI_METHOD_SUPPORT").sum()),
            "nominal_directionally_concordant": int(selected["robustness_status"].eq("NOMINAL_DIRECTIONALLY_CONCORDANT").sum()),
            "nominal_direction_discordance": int(selected["robustness_status"].eq("NOMINAL_IVW_ONLY_DIRECTION_DISCORDANCE").sum()),
            "selected_taxa_in_identical_iv_set_groups": int(selected["non_independent_identical_iv_set"].eq("yes").sum()),
        },
        "gutmgene_raw_encoding": raw_encoding,
        "evidence_boundary": {
            "strict": "Only same-rank normalized exact matches are treated as strict coverage.",
            "permissive": "Genus-to-species matches are an exploratory upper bound and are not merged with strict evidence.",
            "causality": "gutMGene associative mode is retained verbatim; correlative edges are not interpreted as microbial production or causal mediation.",
        },
    }
    (qa_dir / "stage2_nominal_selection_gutmgene_metrics.json").write_text(
        json.dumps(metrics, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    print(
        "PASS_NOMINAL_SELECTION_WITH_LIMITED_STRICT_GUTMGENE_COVERAGE "
        f"selected={len(selected)} strict_taxa={metrics['counts']['strict_covered_selected_taxa']} "
        f"strict_edges={len(strict_selected)} permissive_edges={len(permissive)}"
    )


if __name__ == "__main__":
    main()
