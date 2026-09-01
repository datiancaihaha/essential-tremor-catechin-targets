import argparse
import json
import platform
import sys
import time
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import patsy
import tensorflow as tf
import tensorflow_probability as tfp
from sccoda.util import comp_ana


METADATA_COLUMNS = {
    "donor_id",
    "disease",
    "sex",
    "age",
    "batch_profile",
    "condition",
    "age_c",
}


def select_references(counts, cell_types, number=3):
    percent_zero = np.mean(counts == 0, axis=0)
    relative = counts / counts.sum(axis=1, keepdims=True)
    dispersion = np.var(relative, axis=0) / np.mean(relative, axis=0)
    eligible = np.where(percent_zero < 0.05)[0]
    if len(eligible) < number:
        raise ValueError(
            f"Only {len(eligible)} cell types are present in more than 95% of samples"
        )
    ordered = eligible[np.argsort(dispersion[eligible])]
    table = pd.DataFrame(
        {
            "cell_type": cell_types,
            "fraction_zero": percent_zero,
            "relative_abundance_dispersion": dispersion,
            "eligible_reference": percent_zero < 0.05,
        }
    ).sort_values("relative_abundance_dispersion", kind="stable")
    return [cell_types[index] for index in ordered[:number]], table


def prepare_dataset(path):
    frame = pd.read_csv(path, dtype={"donor_id": str})
    cell_types = [column for column in frame.columns if column not in METADATA_COLUMNS]
    if len(cell_types) != 14:
        raise ValueError(f"Expected 14 cell types; found {len(cell_types)}")
    counts = frame[cell_types].to_numpy(dtype=np.float64)
    if np.any(counts < 0) or not np.allclose(counts, np.round(counts)):
        raise ValueError("scCODA input must contain non-negative integer counts")
    if np.any(counts.sum(axis=1) == 0):
        raise ValueError("Each donor must have at least one observed nucleus")
    observation = frame[
        ["donor_id", "condition", "age_c", "sex", "batch_profile"]
    ].copy()
    observation["condition"] = pd.Categorical(
        observation["condition"], categories=["control", "case"]
    )
    observation["sex"] = pd.Categorical(
        observation["sex"], categories=["female", "male"]
    )
    observation["batch_profile"] = pd.Categorical(observation["batch_profile"])
    formula = (
        "C(condition, Treatment('control')) + age_c + "
        "C(sex, Treatment('female')) + C(batch_profile)"
    )
    design = patsy.dmatrix(formula, observation)
    rank = np.linalg.matrix_rank(np.asarray(design))
    if rank != design.shape[1]:
        raise ValueError(
            f"scCODA design matrix is not full rank: rank={rank}, columns={design.shape[1]}"
        )
    data = ad.AnnData(
        X=counts,
        obs=observation.set_index("donor_id"),
        var=pd.DataFrame(index=cell_types),
    )
    return data, formula, rank, design.shape[1]


def flatten_effect_table(result, analysis_set, requested_reference, actual_reference):
    _, effect = result.summary_prepare(est_fdr=0.05, hdi_prob=0.95)
    effect = effect.reset_index()
    effect["analysis_set"] = analysis_set
    effect["requested_reference"] = requested_reference
    effect["reference_cell_type"] = actual_reference
    effect["credible_effect_at_sccoda_estimated_fdr_0_05"] = (
        effect["Final Parameter"] != 0
    )
    effect["sccoda_inclusion_probability_threshold"] = result.model_specs.get(
        "threshold_prob", np.nan
    )
    return effect


def run_models(input_path, analysis_set, out_dir, seed_offset, sampler):
    data, formula, design_rank, design_columns = prepare_dataset(input_path)
    cell_types = list(data.var_names)
    references, reference_table = select_references(
        np.asarray(data.X), cell_types, number=3
    )
    reference_table["analysis_set"] = analysis_set
    reference_table.to_csv(
        out_dir / f"{analysis_set}_reference_selection.csv", index=False
    )

    requested_references = ["automatic", references[1], references[2]]
    effect_tables = []
    model_audit = []
    for model_index, requested_reference in enumerate(requested_references):
        seed = 20260828 + seed_offset + model_index
        np.random.seed(seed)
        tf.random.set_seed(seed)
        model = comp_ana.CompositionalAnalysis(
            data,
            formula=formula,
            reference_cell_type=requested_reference,
            automatic_reference_absence_threshold=0.05,
        )
        actual_reference = model.cell_types[model.reference_cell_type]
        started = time.time()
        sample_method = model.sample_hmc_da if sampler == "hmc_da" else model.sample_hmc
        result = sample_method(
            num_results=20000,
            num_burnin=5000,
            num_adapt_steps=4000,
            num_leapfrog_steps=10,
            step_size=0.01,
            verbose=False,
        )
        duration = time.time() - started
        effect = flatten_effect_table(
            result,
            analysis_set,
            requested_reference,
            actual_reference,
        )
        effect_tables.append(effect)
        condition_effects = effect[
            effect["Covariate"].astype(str).str.contains(
                "condition", case=False, regex=False
            )
        ]
        model_audit.append(
            {
                "analysis_set": analysis_set,
                "sampler": sampler,
                "requested_reference": requested_reference,
                "reference_cell_type": actual_reference,
                "samples": data.n_obs,
                "cell_types": data.n_vars,
                "design_rank": design_rank,
                "design_columns": design_columns,
                "random_seed": seed,
                "num_results": 20000,
                "num_burnin": 5000,
                "num_adapt_steps": 4000,
                "num_leapfrog_steps": 10,
                "duration_seconds": duration,
                "acceptance_rate": result.sampling_stats.get("acc_rate", np.nan),
                "inclusion_probability_threshold": result.model_specs.get(
                    "threshold_prob", np.nan
                ),
                "credible_condition_effects": int(
                    condition_effects[
                        "credible_effect_at_sccoda_estimated_fdr_0_05"
                    ].sum()
                ),
            }
        )
        print(
            "SCCODA_MODEL_COMPLETE "
            f"analysis_set={analysis_set} "
            f"sampler={sampler} "
            f"reference={actual_reference} "
            f"acceptance_rate={result.sampling_stats.get('acc_rate', np.nan):.4f} "
            f"duration_seconds={duration:.1f}",
            flush=True,
        )

    effects = pd.concat(effect_tables, ignore_index=True)
    audit = pd.DataFrame(model_audit)
    effects.to_csv(out_dir / f"{analysis_set}_all_effects.csv", index=False)
    condition_effects = effects[
        effects["Covariate"].astype(str).str.contains(
            "condition", case=False, regex=False
        )
    ].copy()
    condition_effects.to_csv(
        out_dir / f"{analysis_set}_condition_effects.csv", index=False
    )
    audit.to_csv(out_dir / f"{analysis_set}_model_audit.csv", index=False)
    return condition_effects, audit


def summarize_reference_stability(condition_effects):
    summaries = []
    for (analysis_set, cell_type), values in condition_effects.groupby(
        ["analysis_set", "Cell Type"], sort=True
    ):
        nonreference = values[values["Cell Type"] != values["reference_cell_type"]]
        log_fold = nonreference["log2-fold change"].to_numpy(dtype=float)
        summaries.append(
            {
                "analysis_set": analysis_set,
                "cell_type": cell_type,
                "evaluable_reference_models": len(nonreference),
                "credible_reference_models": int(
                    nonreference[
                        "credible_effect_at_sccoda_estimated_fdr_0_05"
                    ].sum()
                ),
                "same_direction_across_evaluable_references": bool(
                    len(log_fold) > 0
                    and (
                        np.all(log_fold >= 0)
                        or np.all(log_fold <= 0)
                    )
                ),
                "minimum_log2_fold_change": (
                    float(np.min(log_fold)) if len(log_fold) else np.nan
                ),
                "maximum_log2_fold_change": (
                    float(np.max(log_fold)) if len(log_fold) else np.nan
                ),
            }
        )
    return pd.DataFrame(summaries)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--full-input", required=True)
    parser.add_argument("--matched-input", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument(
        "--sampler", choices=["hmc", "hmc_da"], default="hmc"
    )
    args = parser.parse_args()

    out_dir = Path(args.out)
    out_dir.mkdir(parents=True, exist_ok=True)
    full_effects, full_audit = run_models(
        args.full_input, "all_donors", out_dir, seed_offset=0, sampler=args.sampler
    )
    matched_effects, matched_audit = run_models(
        args.matched_input,
        "matched_donors",
        out_dir,
        seed_offset=100,
        sampler=args.sampler,
    )
    all_effects = pd.concat([full_effects, matched_effects], ignore_index=True)
    stability = summarize_reference_stability(all_effects)
    stability.to_csv(out_dir / "sccoda_reference_stability.csv", index=False)

    audit = pd.concat([full_audit, matched_audit], ignore_index=True)
    package_versions = {
        "python": platform.python_version(),
        "numpy": np.__version__,
        "pandas": pd.__version__,
        "anndata": ad.__version__,
        "tensorflow": tf.__version__,
        "tensorflow_probability": tfp.__version__,
        "sccoda": "0.1.9",
        "sampler": args.sampler,
        "models_completed": int(len(audit)),
        "all_models_completed": bool(len(audit) == 6),
    }
    with open(out_dir / "software_and_completion.json", "w", encoding="utf-8") as handle:
        json.dump(package_versions, handle, ensure_ascii=False, indent=2)

    summary_lines = [
        "# SCP3177 scCODA组成敏感性分析",
        "",
        "## 方法",
        "",
        "使用供体层14类细胞原始计数，模型包含ET状态、年龄、性别和测序批次组合。",
        "全供体和年龄/性别匹配供体分别采用自动参考细胞类型及两个低离散度备选参考细胞类型。",
        "scCODA结果按其贝叶斯后验可信效应和估计FDR报告；这些结果仅作为组成敏感性证据，不用于新增、删除或裁剪propeller按nominal P < 0.05确定的候选。",
        "",
        "## 完成度",
        "",
        f"完成{len(audit)}个预设scCODA模型。",
        "",
        "## 方法依据",
        "",
        "- Büttner M, et al. Nature Communications. 2021;12:6876. doi:10.1038/s41467-021-27150-6.",
    ]
    (out_dir / "SCCODA_RESULTS_zh.md").write_text(
        "\n".join(summary_lines) + "\n", encoding="utf-8"
    )


if __name__ == "__main__":
    main()
