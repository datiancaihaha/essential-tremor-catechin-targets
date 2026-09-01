suppressPackageStartupMessages({
  library(data.table)
  library(jsonlite)
  library(TwoSampleMR)
})

args <- commandArgs(trailingOnly = TRUE)
project_root <- if (length(args) >= 1) normalizePath(args[[1]], winslash = "/", mustWork = TRUE) else normalizePath(".", winslash = "/", mustWork = TRUE)
input_path <- file.path(project_root, "outputs", "stage1_clumped_independent_ivs.csv")
output_dir <- file.path(project_root, "outputs")
qa_dir <- file.path(project_root, "qa")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(qa_dir, recursive = TRUE, showWarnings = FALSE)

analysis_seed <- 20260821L
parameters <- TwoSampleMR::default_parameters()
parameters$nboot <- 1000L
parameters$phi <- 1

stop_with <- function(message) stop(message, call. = FALSE)
finite_required <- function(x, label) {
  if (any(!is.finite(x))) stop_with(paste0("Non-finite values in ", label))
}

dat_all <- data.table::fread(input_path, na.strings = c("", "NA", "NaN"))
required_columns <- c(
  "bac", "analysis_rsid", "beta_exposure", "se_exposure", "beta_outcome_aligned",
  "se_outcome", "f_statistic", "mr_usable"
)
missing_columns <- setdiff(required_columns, names(dat_all))
if (length(missing_columns) > 0) stop_with(paste("Missing required columns:", paste(missing_columns, collapse = ", ")))

dat <- dat_all[mr_usable == "yes"]
if (nrow(dat) == 0) stop_with("No mr_usable=yes rows")
finite_required(dat$beta_exposure, "beta_exposure")
finite_required(dat$se_exposure, "se_exposure")
finite_required(dat$beta_outcome_aligned, "beta_outcome_aligned")
finite_required(dat$se_outcome, "se_outcome")
finite_required(dat$f_statistic, "f_statistic")
if (any(dat$se_exposure <= 0 | dat$se_outcome <= 0)) stop_with("Non-positive standard errors detected")
if (any(dat$f_statistic < 10)) stop_with("Weak instrument with F<10 detected in MR input")
if (anyDuplicated(dat[, .(bac, analysis_rsid)])) stop_with("Duplicate taxon-rsID rows detected")

taxa <- sort(unique(dat$bac))
taxon_counts <- dat[, .N, by = bac]
if (length(taxa) != 211L) stop_with(paste0("Expected 211 taxa, found ", length(taxa)))
if (min(taxon_counts$N) < 3L) stop_with("At least one taxon has fewer than 3 usable IVs")

result_row <- function(taxon, method, x) {
  ok <- is.list(x) && all(c("b", "se", "pval") %in% names(x)) &&
    is.finite(x$b) && is.finite(x$se) && x$se > 0 && is.finite(x$pval)
  if (!ok) {
    return(data.table(
      bac = taxon, method = method, nsnp = length(dat[bac == taxon]$analysis_rsid),
      beta = NA_real_, se = NA_real_, p_value = NA_real_, ci_lower = NA_real_, ci_upper = NA_real_,
      odds_ratio = NA_real_, or_ci_lower = NA_real_, or_ci_upper = NA_real_,
      method_status = "METHOD_FAILED", method_status_zh = "方法运行失败或返回非有限值"
    ))
  }
  ci_lower <- x$b - 1.96 * x$se
  ci_upper <- x$b + 1.96 * x$se
  data.table(
    bac = taxon, method = method, nsnp = as.integer(x$nsnp %||% length(dat[bac == taxon]$analysis_rsid)),
    beta = x$b, se = x$se, p_value = x$pval, ci_lower = ci_lower, ci_upper = ci_upper,
    odds_ratio = exp(x$b), or_ci_lower = exp(ci_lower), or_ci_upper = exp(ci_upper),
    method_status = "OK", method_status_zh = "方法成功运行"
  )
}

`%||%` <- function(x, y) if (is.null(x) || length(x) == 0 || all(is.na(x))) y else x

method_results <- vector("list", length(taxa) * 4L)
sensitivity_results <- vector("list", length(taxa))
loo_results <- vector("list", length(taxa))
position <- 1L

for (i in seq_along(taxa)) {
  taxon <- taxa[[i]]
  x <- dat[bac == taxon]
  bx <- x$beta_exposure
  by <- x$beta_outcome_aligned
  sx <- x$se_exposure
  sy <- x$se_outcome

  ivw <- tryCatch(TwoSampleMR::mr_ivw(bx, by, sx, sy, parameters), error = function(e) list(error = conditionMessage(e)))
  set.seed(analysis_seed + i)
  weighted_median <- tryCatch(TwoSampleMR::mr_weighted_median(bx, by, sx, sy, parameters), error = function(e) list(error = conditionMessage(e)))
  egger <- tryCatch(TwoSampleMR::mr_egger_regression(bx, by, sx, sy, parameters), error = function(e) list(error = conditionMessage(e)))
  set.seed(analysis_seed + 100000L + i)
  weighted_mode <- tryCatch(TwoSampleMR::mr_weighted_mode(bx, by, sx, sy, parameters), error = function(e) list(error = conditionMessage(e)))

  method_results[[position]] <- result_row(taxon, "Inverse variance weighted", ivw)
  method_results[[position + 1L]] <- result_row(taxon, "Weighted median", weighted_median)
  method_results[[position + 2L]] <- result_row(taxon, "MR Egger", egger)
  method_results[[position + 3L]] <- result_row(taxon, "Weighted mode", weighted_mode)
  position <- position + 4L

  ivw_i2 <- if (is.finite(ivw$Q %||% NA_real_) && ivw$Q > 0) max(0, (ivw$Q - ivw$Q_df) / ivw$Q) else NA_real_
  egger_i2 <- if (is.finite(egger$Q %||% NA_real_) && egger$Q > 0) max(0, (egger$Q - egger$Q_df) / egger$Q) else NA_real_
  sensitivity_results[[i]] <- data.table(
    bac = taxon,
    nsnp = nrow(x),
    f_min = min(x$f_statistic),
    f_median = median(x$f_statistic),
    f_mean = mean(x$f_statistic),
    ivw_q = ivw$Q %||% NA_real_,
    ivw_q_df = ivw$Q_df %||% NA_real_,
    ivw_q_p = ivw$Q_pval %||% NA_real_,
    ivw_i2 = ivw_i2,
    egger_q = egger$Q %||% NA_real_,
    egger_q_df = egger$Q_df %||% NA_real_,
    egger_q_p = egger$Q_pval %||% NA_real_,
    egger_i2 = egger_i2,
    egger_intercept = egger$b_i %||% NA_real_,
    egger_intercept_se = egger$se_i %||% NA_real_,
    egger_intercept_p = egger$pval_i %||% NA_real_
  )

  full_beta <- ivw$b %||% NA_real_
  taxon_loo <- rbindlist(lapply(seq_len(nrow(x)), function(j) {
    loo <- tryCatch(
      TwoSampleMR::mr_ivw(bx[-j], by[-j], sx[-j], sy[-j], parameters),
      error = function(e) list(error = conditionMessage(e))
    )
    ok <- is.list(loo) && all(c("b", "se", "pval") %in% names(loo)) &&
      is.finite(loo$b) && is.finite(loo$se) && loo$se > 0 && is.finite(loo$pval)
    if (!ok) {
      return(data.table(
        bac = taxon, omitted_rsid = x$analysis_rsid[[j]], remaining_nsnp = nrow(x) - 1L,
        beta = NA_real_, se = NA_real_, p_value = NA_real_, ci_lower = NA_real_, ci_upper = NA_real_,
        odds_ratio = NA_real_, or_ci_lower = NA_real_, or_ci_upper = NA_real_, delta_beta_vs_full = NA_real_,
        method_status = "METHOD_FAILED", method_status_zh = "留一法IVW运行失败或返回非有限值"
      ))
    }
    ci_lower <- loo$b - 1.96 * loo$se
    ci_upper <- loo$b + 1.96 * loo$se
    data.table(
      bac = taxon, omitted_rsid = x$analysis_rsid[[j]], remaining_nsnp = nrow(x) - 1L,
      beta = loo$b, se = loo$se, p_value = loo$pval, ci_lower = ci_lower, ci_upper = ci_upper,
      odds_ratio = exp(loo$b), or_ci_lower = exp(ci_lower), or_ci_upper = exp(ci_upper),
      delta_beta_vs_full = loo$b - full_beta,
      method_status = "OK", method_status_zh = "留一法IVW成功运行"
    )
  }))
  loo_results[[i]] <- taxon_loo
}

methods <- rbindlist(method_results, use.names = TRUE, fill = TRUE)
methods[, method_order := match(method, c("Inverse variance weighted", "Weighted median", "MR Egger", "Weighted mode"))]
setorder(methods, bac, method_order)
methods[, method_order := NULL]

ivw_primary <- methods[method == "Inverse variance weighted"]
if (nrow(ivw_primary) != 211L || any(ivw_primary$method_status != "OK")) stop_with("IVW primary analysis did not complete for all 211 taxa")
ivw_primary[, fdr_bh := p.adjust(p_value, method = "BH")]
ivw_primary[, evidence_status := fifelse(
  fdr_bh < 0.05,
  "FDR_SIGNIFICANT",
  fifelse(p_value < 0.05, "NOMINAL_SUGGESTIVE", "NO_PRIMARY_EVIDENCE")
)]
ivw_primary[, evidence_status_zh := fifelse(
  evidence_status == "FDR_SIGNIFICANT",
  "FDR显著（211个IVW主检验经BH校正后q<0.05）",
  fifelse(
    evidence_status == "NOMINAL_SUGGESTIVE",
    "名义提示性（IVW P<0.05，但未通过BH-FDR）",
    "主分析未见统计学证据（IVW P>=0.05）"
  )
)]
setorder(ivw_primary, p_value, bac)

sensitivity <- rbindlist(sensitivity_results, use.names = TRUE, fill = TRUE)
loo <- rbindlist(loo_results, use.names = TRUE, fill = TRUE)
loo_summary <- loo[, {
  ok <- method_status == "OK" & is.finite(beta)
  full <- ivw_primary[bac == .BY$bac]$beta[[1]]
  if (!any(ok)) {
    list(
      loo_runs = .N, loo_successful = 0L, loo_beta_min = NA_real_, loo_beta_max = NA_real_,
      max_abs_delta_beta = NA_real_, most_influential_omitted_rsid = NA_character_,
      any_direction_flip = NA, all_loo_direction_concordant = NA
    )
  } else {
    z <- .SD[ok]
    k <- which.max(abs(z$delta_beta_vs_full))
    flips <- if (full == 0) rep(NA, nrow(z)) else sign(z$beta) != sign(full)
    list(
      loo_runs = .N, loo_successful = nrow(z), loo_beta_min = min(z$beta), loo_beta_max = max(z$beta),
      max_abs_delta_beta = max(abs(z$delta_beta_vs_full)), most_influential_omitted_rsid = z$omitted_rsid[[k]],
      any_direction_flip = if (all(is.na(flips))) NA else any(flips, na.rm = TRUE),
      all_loo_direction_concordant = if (all(is.na(flips))) NA else !any(flips, na.rm = TRUE)
    )
  }
}, by = bac]

method_wide <- dcast(methods, bac ~ method, value.var = c("beta", "p_value", "method_status"))
setnames(
  method_wide,
  old = names(method_wide),
  new = gsub("[^A-Za-z0-9]+", "_", tolower(names(method_wide)))
)
decision <- merge(ivw_primary, method_wide, by = "bac", all.x = TRUE, sort = FALSE)
decision <- merge(decision, sensitivity, by = "bac", all.x = TRUE, sort = FALSE)
decision <- merge(decision, loo_summary, by = "bac", all.x = TRUE, sort = FALSE)
setorder(decision, p_value, bac)

data.table::fwrite(methods, file.path(output_dir, "stage2_mr_results_all_methods.csv"), na = "")
data.table::fwrite(ivw_primary, file.path(output_dir, "stage2_ivw_primary_results.csv"), na = "")
data.table::fwrite(sensitivity, file.path(output_dir, "stage2_sensitivity_summary.csv"), na = "")
data.table::fwrite(loo, file.path(output_dir, "stage2_leave_one_out.csv"), na = "")
data.table::fwrite(loo_summary, file.path(output_dir, "stage2_leave_one_out_summary.csv"), na = "")
data.table::fwrite(decision, file.path(output_dir, "stage2_taxon_core_decision_summary.csv"), na = "")

metrics <- list(
  status = "PASS_CORE_MR_EXECUTION",
  status_zh = "通过：211个taxa的核心MR及敏感性分析已执行",
  generated_at_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
  analysis_seed = analysis_seed,
  input = list(path = input_path, selected_iv_rows = nrow(dat), taxa = length(taxa), min_iv_per_taxon = min(taxon_counts$N), max_iv_per_taxon = max(taxon_counts$N)),
  specification = list(
    primary_method = "TwoSampleMR::mr_ivw (multiplicative random effects with under-dispersion correction)",
    sensitivity_methods = c("TwoSampleMR::mr_weighted_median", "TwoSampleMR::mr_egger_regression", "TwoSampleMR::mr_weighted_mode"),
    bootstrap_repetitions = parameters$nboot,
    weighted_mode_phi = parameters$phi,
    multiplicity_family = "211 IVW primary p-values",
    multiplicity_method = "Benjamini-Hochberg FDR",
    fdr_threshold = 0.05,
    nominal_threshold = 0.05
  ),
  results = list(
    method_rows = nrow(methods),
    method_failed_rows = sum(methods$method_status != "OK"),
    ivw_rows = nrow(ivw_primary),
    fdr_significant_taxa = sum(ivw_primary$fdr_bh < 0.05),
    nominal_ivw_taxa = sum(ivw_primary$p_value < 0.05),
    nominal_only_taxa = sum(ivw_primary$p_value < 0.05 & ivw_primary$fdr_bh >= 0.05),
    loo_rows = nrow(loo),
    loo_failed_rows = sum(loo$method_status != "OK"),
    taxa_with_any_loo_direction_flip = sum(loo_summary$any_direction_flip %in% TRUE)
  ),
  software = list(
    R = R.version.string,
    TwoSampleMR = as.character(utils::packageVersion("TwoSampleMR")),
    data_table = as.character(utils::packageVersion("data.table")),
    jsonlite = as.character(utils::packageVersion("jsonlite"))
  )
)
jsonlite::write_json(metrics, file.path(qa_dir, "stage2_core_mr_metrics.json"), pretty = TRUE, auto_unbox = TRUE, na = "null")

cat(sprintf("PASS_CORE_MR_EXECUTION taxa=%d iv_rows=%d fdr=%d nominal=%d method_failures=%d loo_failures=%d\n",
  length(taxa), nrow(dat), metrics$results$fdr_significant_taxa, metrics$results$nominal_ivw_taxa,
  metrics$results$method_failed_rows, metrics$results$loo_failed_rows
))
