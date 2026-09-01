suppressPackageStartupMessages({
  library(data.table)
  library(jsonlite)
})

args <- commandArgs(trailingOnly = TRUE)
project_root <- if (length(args) >= 1) normalizePath(args[[1]], winslash = "/", mustWork = TRUE) else normalizePath(".", winslash = "/", mustWork = TRUE)
.libPaths(c(file.path(project_root, "tools", "R_library"), .libPaths()))
suppressPackageStartupMessages({
  library(TwoSampleMR)
  library(mr.raps)
})

input_path <- file.path(project_root, "outputs", "stage1_clumped_independent_ivs.csv")
selected_path <- file.path(project_root, "outputs", "stage2_nominal_selected_taxa.csv")
output_dir <- file.path(project_root, "outputs")
qa_dir <- file.path(project_root, "qa")
report_path <- file.path(project_root, "STAGE5_MR_STEIGER_RAPS_DECISION_zh.md")

priority_taxa <- c(
  "genus.Faecalibacterium.id.2057",
  "genus.Flavonifractor.id.2059",
  "genus.RuminococcaceaeUCG011.id.11368",
  "genus.Methanobrevibacter.id.123"
)
prevalence_grid <- c(0.0032, 0.0133, 0.0579)
prevalence_labels <- c("all_age_meta_analysis_0.32pct", "all_age_meta_analysis_1.33pct", "age_65plus_meta_analysis_5.79pct")
ncase_et <- 16480
ncontrol_et <- 1936173

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(qa_dir, recursive = TRUE, showWarnings = FALSE)

stop_with <- function(message) stop(message, call. = FALSE)

dat <- data.table::fread(input_path, na.strings = c("", "NA", "NaN"))
selected <- data.table::fread(selected_path, na.strings = c("", "NA", "NaN"))
dat <- dat[bac %in% priority_taxa & mr_usable == "yes"]
if (!setequal(unique(dat$bac), priority_taxa)) stop_with("One or more priority taxa lack usable IVs")
if (any(!is.finite(dat$beta_exposure)) || any(!is.finite(dat$se_exposure)) ||
    any(!is.finite(dat$beta_outcome_aligned)) || any(!is.finite(dat$se_outcome)) ||
    any(!is.finite(dat$outcome_eaf_median)) || any(!is.finite(dat$n_exposure))) {
  stop_with("Non-finite values in MR-Steiger/MR-RAPS input")
}
if (any(dat$outcome_eaf_median <= 0 | dat$outcome_eaf_median >= 1)) stop_with("Invalid outcome allele frequency")

capture_warnings <- function(expression) {
  messages <- character()
  value <- withCallingHandlers(
    expression,
    warning = function(warning) {
      messages <<- c(messages, conditionMessage(warning))
      invokeRestart("muffleWarning")
    }
  )
  list(value = value, warnings = unique(messages))
}

steiger_rows <- list()
steiger_snp_rows <- list()
steiger_index <- 1L
snp_index <- 1L
for (taxon in priority_taxa) {
  x <- dat[bac == taxon]
  taxon_name <- selected[bac == taxon]$taxon_name[[1]]
  r_exposure <- TwoSampleMR::get_r_from_bsen(x$beta_exposure, x$se_exposure, x$n_exposure)
  if (any(!is.finite(r_exposure))) stop_with(paste("Non-finite exposure r for", taxon))
  for (j in seq_along(prevalence_grid)) {
    prevalence <- prevalence_grid[[j]]
    r_outcome <- TwoSampleMR::get_r_from_lor(
      lor = x$beta_outcome_aligned,
      af = x$outcome_eaf_median,
      ncase = rep(ncase_et, nrow(x)),
      ncontrol = rep(ncontrol_et, nrow(x)),
      prevalence = rep(prevalence, nrow(x)),
      model = "logit",
      correction = FALSE
    )
    if (any(!is.finite(r_outcome))) stop_with(paste("Non-finite outcome r for", taxon, prevalence))
    direction_input <- data.frame(
      id.exposure = taxon,
      id.outcome = "deCODE_ET_G250",
      exposure = taxon_name,
      outcome = "Essential tremor",
      pval.exposure = x$p_exposure,
      pval.outcome = x$p_outcome,
      samplesize.exposure = x$n_exposure,
      samplesize.outcome = ncase_et + ncontrol_et,
      r.exposure = r_exposure,
      r.outcome = r_outcome,
      stringsAsFactors = FALSE
    )
    result <- TwoSampleMR::directionality_test(direction_input)
    if (nrow(result) != 1L) stop_with(paste("Unexpected Steiger result count for", taxon))
    steiger_rows[[steiger_index]] <- data.table(
      bac = taxon,
      taxon_name = taxon_name,
      nsnp = nrow(x),
      et_population_prevalence = prevalence,
      prevalence_scenario = prevalence_labels[[j]],
      snp_r2_exposure = result$snp_r2.exposure,
      snp_r2_outcome = result$snp_r2.outcome,
      r2_exposure_to_outcome_ratio = result$snp_r2.exposure / result$snp_r2.outcome,
      correct_causal_direction = result$correct_causal_direction,
      steiger_p_value = result$steiger_pval,
      interpretation_zh = ifelse(
        result$correct_causal_direction,
        "工具变量解释的菌群方差大于ET责任度方差，方向与菌群到ET一致",
        "未支持预设的菌群到ET方向"
      )
    )
    steiger_index <- steiger_index + 1L
    for (k in seq_len(nrow(x))) {
      steiger_snp_rows[[snp_index]] <- data.table(
        bac = taxon,
        taxon_name = taxon_name,
        analysis_rsid = x$analysis_rsid[[k]],
        et_population_prevalence = prevalence,
        prevalence_scenario = prevalence_labels[[j]],
        r_exposure = r_exposure[[k]],
        r2_exposure = r_exposure[[k]]^2,
        r_outcome = r_outcome[[k]],
        r2_outcome = r_outcome[[k]]^2,
        snp_supports_exposure_to_outcome = r_exposure[[k]]^2 > r_outcome[[k]]^2
      )
      snp_index <- snp_index + 1L
    }
  }
}
steiger <- data.table::rbindlist(steiger_rows)
steiger_snp <- data.table::rbindlist(steiger_snp_rows)

raps_rows <- list()
for (i in seq_along(priority_taxa)) {
  taxon <- priority_taxa[[i]]
  x <- dat[bac == taxon]
  taxon_name <- selected[bac == taxon]$taxon_name[[1]]
  ivw <- selected[bac == taxon]
  overdispersed <- capture_warnings(
    mr.raps::mr.raps.overdispersed.robust(
      x$beta_exposure,
      x$beta_outcome_aligned,
      x$se_exposure,
      x$se_outcome,
      loss.function = "huber",
      initialization = "l2",
      suppress.warning = FALSE,
      diagnosis = FALSE
    )
  )
  simple <- capture_warnings(
    mr.raps::mr.raps.simple.robust(
      x$beta_exposure,
      x$beta_outcome_aligned,
      x$se_exposure,
      x$se_outcome,
      loss.function = "huber",
      diagnosis = FALSE
    )
  )
  fit <- overdispersed$value
  simple_fit <- simple$value
  raps_rows[[i]] <- data.table(
    bac = taxon,
    taxon_name = taxon_name,
    nsnp = nrow(x),
    f_min = min(x$f_statistic),
    f_median = median(x$f_statistic),
    ivw_beta = ivw$beta[[1]],
    ivw_p_value = ivw$p_value[[1]],
    raps_model = "overdispersed robust adjusted profile score; Huber loss; sandwich SE",
    raps_beta = fit$beta.hat,
    raps_se = fit$beta.se,
    raps_p_value = fit$beta.p.value,
    raps_odds_ratio = exp(fit$beta.hat),
    raps_or_ci_lower = exp(fit$beta.hat - 1.96 * fit$beta.se),
    raps_or_ci_upper = exp(fit$beta.hat + 1.96 * fit$beta.se),
    raps_tau2 = fit$tau2.hat,
    raps_tau2_se = fit$tau2.se,
    raps_warning = paste(overdispersed$warnings, collapse = " | "),
    simple_robust_beta = simple_fit$beta.hat,
    simple_robust_se = simple_fit$beta.se,
    simple_robust_p_value = simple_fit$beta.p.value,
    simple_robust_odds_ratio = exp(simple_fit$beta.hat),
    simple_robust_or_ci_lower = exp(simple_fit$beta.hat - 1.96 * simple_fit$beta.se),
    simple_robust_or_ci_upper = exp(simple_fit$beta.hat + 1.96 * simple_fit$beta.se),
    simple_robust_warning = paste(simple$warnings, collapse = " | "),
    direction_concordant_with_ivw = sign(fit$beta.hat) == sign(ivw$beta[[1]]),
    analysis_role_zh = "稳健性敏感性分析；不改变nominal P<0.05入选集合"
  )
}
raps <- data.table::rbindlist(raps_rows)

data.table::fwrite(steiger, file.path(output_dir, "stage5_mr_steiger_prevalence_sensitivity.csv"), na = "")
data.table::fwrite(steiger_snp, file.path(output_dir, "stage5_mr_steiger_snp_r2_audit.csv"), na = "")
data.table::fwrite(raps, file.path(output_dir, "stage5_mr_raps_priority_taxa.csv"), na = "")

all_direction_supported <- all(steiger$correct_causal_direction)
all_raps_concordant <- all(raps$direction_concordant_with_ivw)
checks <- list(
  four_priority_taxa_present = setequal(unique(dat$bac), priority_taxa),
  steiger_three_prevalence_scenarios_per_taxon = nrow(steiger) == length(priority_taxa) * length(prevalence_grid),
  steiger_all_values_finite = all(is.finite(steiger$snp_r2_exposure)) && all(is.finite(steiger$snp_r2_outcome)) && all(is.finite(steiger$steiger_p_value)),
  steiger_snp_audit_complete = nrow(steiger_snp) == sum(dat[, .N, by = bac]$N) * length(prevalence_grid),
  mr_raps_four_results = nrow(raps) == length(priority_taxa),
  mr_raps_all_values_finite = all(is.finite(raps$raps_beta)) && all(is.finite(raps$raps_se)) && all(is.finite(raps$raps_p_value)),
  mr_raps_simple_sensitivity_all_values_finite = all(is.finite(raps$simple_robust_beta)) && all(is.finite(raps$simple_robust_se)) && all(is.finite(raps$simple_robust_p_value)),
  nominal_selection_rule_unchanged = TRUE
)
qa_status <- if (all(unlist(checks))) "PASS_STAGE5_MR_STEIGER_RAPS_QA" else "FAIL_STAGE5_MR_STEIGER_RAPS_QA"
qa_status_zh <- if (startsWith(qa_status, "PASS")) "通过：MR-Steiger患病率敏感性分析及4个重点菌MR-RAPS均完整运行。" else "失败：MR-Steiger或MR-RAPS完整性检查未通过。"
Encoding(qa_status_zh) <- "UTF-8"
qa <- list(
  generated_at_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
  status = qa_status,
  status_zh = qa_status_zh,
  checks = checks,
  results_not_qa_gates = list(
    all_prevalence_scenarios_support_exposure_to_outcome = all_direction_supported,
    all_mr_raps_directions_concordant_with_ivw = all_raps_concordant
  ),
  parameters = list(
    et_cases = ncase_et,
    et_controls = ncontrol_et,
    et_population_prevalence_sensitivity = prevalence_grid,
    steiger_exposure_r = "TwoSampleMR::get_r_from_bsen",
    steiger_binary_outcome_r = "TwoSampleMR::get_r_from_lor",
    mr_raps_primary = "mr.raps::mr.raps.overdispersed.robust with Huber loss",
    mr_raps_secondary = "mr.raps::mr.raps.simple.robust with Huber loss"
  ),
  software = list(
    R = R.version.string,
    TwoSampleMR = as.character(utils::packageVersion("TwoSampleMR")),
    mr_raps = as.character(utils::packageVersion("mr.raps"))
  )
)
jsonlite::write_json(qa, file.path(qa_dir, "stage5_mr_steiger_raps_qa.json"), pretty = TRUE, auto_unbox = TRUE, digits = 16)
if (!all(unlist(checks))) stop_with(jsonlite::toJSON(qa, auto_unbox = TRUE, pretty = TRUE))

format_bool_zh <- function(value) ifelse(value, "支持预设方向", "不支持预设方向")
report <- c(
  "# ET项目 Stage 5：MR-Steiger与MR-RAPS稳健性分析",
  "",
  paste0("生成时间（UTC）：", format(Sys.time(), tz = "UTC", usetz = TRUE)),
  "",
  "## 判定报告",
  "",
  "`PASS_STAGE5_MR_STEIGER_RAPS_QA`（通过：MR-Steiger患病率敏感性分析及4个重点菌MR-RAPS均完整运行）。",
  "",
  "这些分析是对既有nominal结果的方向性与弱工具变量稳健性核查，不以MR-RAPS或Steiger P值重新筛选、加入或剔除菌群；原IVW nominal P<0.05入选规则不变。",
  "",
  "## MR-Steiger",
  "",
  "ET为二分类结局，责任度相关系数依赖人群患病率。主情景使用全年龄0.32%，并以1.33%和≥65岁5.79%作敏感性范围。",
  "",
  "|菌群|0.32%|1.33%|5.79%|",
  "|---|---|---|---|"
)
for (taxon in priority_taxa) {
  x <- steiger[bac == taxon][order(et_population_prevalence)]
  report <- c(report, paste0(
    "|", x$taxon_name[[1]], "|",
    paste(format_bool_zh(x$correct_causal_direction), collapse = "|"), "|"
  ))
}
report <- c(
  report,
  "",
  paste0("总体方向结果：", ifelse(all_direction_supported, "4个重点菌在全部患病率情景下均支持菌群→ET方向。", "至少一个菌或患病率情景未支持预设方向，详见CSV。")),
  "",
  "MR-Steiger依赖测量精度和责任度尺度假设，只能降低反向解释疑虑，不能替代反向MR或证明因果方向。",
  "",
  "## MR-RAPS",
  "",
  "|菌群|SNP数|MR-RAPS OR (95% CI)|P|与IVW方向一致|",
  "|---|---:|---:|---:|---|"
)
for (i in seq_len(nrow(raps))) {
  row <- raps[i]
  report <- c(report, sprintf(
    "|%s|%d|%.3f (%.3f–%.3f)|%.4g|%s|",
    row$taxon_name, row$nsnp, row$raps_odds_ratio, row$raps_or_ci_lower,
    row$raps_or_ci_upper, row$raps_p_value,
    ifelse(row$direction_concordant_with_ivw, "是", "否")
  ))
}
report <- c(
  report,
  "",
  paste0("方向一致性：", ifelse(all_raps_concordant, "4个重点菌的MR-RAPS方向均与IVW一致。", "至少一个重点菌方向不一致，详见CSV。")),
  "",
  "当过度离散τ²估计接近0时，同时报告simple robust Huber模型；这属于预设的模型敏感性核查，不选择性挑选P值。",
  "",
  "## 证据边界",
  "",
  "- MR-Steiger和MR-RAPS提高的是方向性与稳健性审计完整度，不把nominal关联升级为已证实因果。",
  "- Faecalibacterium仍是核心结果；患者菌群和SCFA公开队列重新分析尚未完成。",
  "- BH-FDR继续作为补充证据强度信息，不用于入选或剔除。",
  ""
)
writeLines(report, report_path, useBytes = TRUE)
cat(jsonlite::toJSON(qa, auto_unbox = TRUE, pretty = TRUE, digits = 16), "\n")
