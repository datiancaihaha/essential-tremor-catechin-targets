suppressPackageStartupMessages({
  library(data.table)
  library(jsonlite)
  library(parallel)
  library(MRPRESSO)
})

args <- commandArgs(trailingOnly = TRUE)
project_root <- if (length(args) >= 1) normalizePath(args[[1]], winslash = "/", mustWork = TRUE) else normalizePath(".", winslash = "/", mustWork = TRUE)
n_workers <- if (length(args) >= 2) as.integer(args[[2]]) else 10L
if (!is.finite(n_workers) || n_workers < 1L) stop("n_workers must be a positive integer", call. = FALSE)

input_path <- file.path(project_root, "outputs", "stage1_clumped_independent_ivs.csv")
output_dir <- file.path(project_root, "outputs")
qa_dir <- file.path(project_root, "qa")
work_dir <- file.path(project_root, "work", "stage2_mr_presso", "run_v1_nb1000")
result_dir <- file.path(work_dir, "taxon_results")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(qa_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)

nb_distribution <- 1000L
significance_threshold <- 0.05
seed_base <- 20260821L

dat <- data.table::fread(input_path, na.strings = c("", "NA", "NaN"))
dat <- dat[mr_usable == "yes"]
taxa <- sort(unique(dat$bac))
if (length(taxa) != 211L) stop(sprintf("Expected 211 taxa, found %d", length(taxa)), call. = FALSE)

result_path <- function(index) file.path(result_dir, sprintf("taxon_%03d.rds", index))

valid_cached_result <- function(index) {
  path <- result_path(index)
  if (!file.exists(path)) return(FALSE)
  z <- tryCatch(readRDS(path), error = function(e) NULL)
  is.list(z) && identical(z$bac, taxa[[index]]) && identical(as.integer(z$nb_distribution), nb_distribution) &&
    identical(as.integer(z$seed), seed_base + index) && !is.null(z$run_status)
}

run_taxon <- function(index) {
  taxon <- taxa[[index]]
  x <- as.data.frame(dat[bac == taxon])
  seed <- seed_base + index
  path <- result_path(index)
  if (nrow(x) <= 3L) {
    z <- list(
      index = index, bac = taxon, nsnp = nrow(x), nb_distribution = nb_distribution,
      seed = seed, run_status = "NOT_APPLICABLE_NOT_ENOUGH_IVS",
      warnings = character(), error = NULL, result = NULL, rsids = x$analysis_rsid
    )
    saveRDS(z, path)
    return(z$run_status)
  }

  rownames(x) <- x$analysis_rsid
  warning_messages <- character()
  set.seed(seed)
  result <- withCallingHandlers(
    tryCatch(
      MRPRESSO::mr_presso(
        BetaOutcome = "beta_outcome_aligned",
        BetaExposure = "beta_exposure",
        SdOutcome = "se_outcome",
        SdExposure = "se_exposure",
        OUTLIERtest = TRUE,
        DISTORTIONtest = TRUE,
        data = x,
        NbDistribution = nb_distribution,
        SignifThreshold = significance_threshold,
        seed = seed
      ),
      error = function(e) structure(list(message = conditionMessage(e)), class = "mrpresso_error")
    ),
    warning = function(w) {
      warning_messages <<- c(warning_messages, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )

  if (inherits(result, "mrpresso_error")) {
    z <- list(
      index = index, bac = taxon, nsnp = nrow(x), nb_distribution = nb_distribution,
      seed = seed, run_status = "ERROR", warnings = unique(warning_messages),
      error = result$message, result = NULL, rsids = x$analysis_rsid
    )
  } else {
    z <- list(
      index = index, bac = taxon, nsnp = nrow(x), nb_distribution = nb_distribution,
      seed = seed, run_status = "OK", warnings = unique(warning_messages),
      error = NULL, result = result, rsids = x$analysis_rsid
    )
  }
  saveRDS(z, path)
  z$run_status
}

pending <- which(!vapply(seq_along(taxa), valid_cached_result, logical(1)))
if (length(pending) > 0L) {
  workers_to_use <- min(n_workers, length(pending))
  cluster_log <- file.path(work_dir, "parallel_workers.log")
  cl <- parallel::makeCluster(workers_to_use, outfile = cluster_log)
  parallel::clusterEvalQ(cl, {
    suppressPackageStartupMessages({
      library(data.table)
      library(MRPRESSO)
    })
    NULL
  })
  parallel::clusterExport(
    cl,
    varlist = c(
      "taxa", "dat", "seed_base", "nb_distribution", "significance_threshold",
      "result_dir", "result_path", "run_taxon"
    ),
    envir = environment()
  )
  invisible(tryCatch(
    parallel::parLapplyLB(cl, pending, run_taxon),
    finally = parallel::stopCluster(cl)
  ))
}

missing_results <- which(!file.exists(vapply(seq_along(taxa), result_path, character(1))))
if (length(missing_results) > 0L) stop(sprintf("Missing %d taxon result files", length(missing_results)), call. = FALSE)
all_results <- lapply(seq_along(taxa), function(i) readRDS(result_path(i)))

parse_p <- function(value) {
  if (is.null(value) || length(value) == 0L || all(is.na(value))) {
    return(list(text = NA_character_, numeric = NA_real_, is_bound = NA))
  }
  text <- as.character(value[[1]])
  is_bound <- startsWith(text, "<")
  numeric <- suppressWarnings(as.numeric(sub("^<", "", text)))
  list(text = text, numeric = numeric, is_bound = is_bound)
}

main_row <- function(main, analysis_name) {
  if (is.null(main) || !is.data.frame(main) || !"MR Analysis" %in% names(main)) return(NULL)
  z <- main[main[["MR Analysis"]] == analysis_name, , drop = FALSE]
  if (nrow(z) == 0L) NULL else z[1, , drop = FALSE]
}

summary_rows <- vector("list", length(all_results))
outlier_rows <- list()
outlier_position <- 1L

for (i in seq_along(all_results)) {
  z <- all_results[[i]]
  if (z$run_status == "NOT_APPLICABLE_NOT_ENOUGH_IVS") {
    summary_rows[[i]] <- data.table(
      bac = z$bac, nsnp = z$nsnp, nb_distribution = z$nb_distribution, seed = z$seed,
      run_status = z$run_status, run_status_zh = "不适用：MR-PRESSO要求至少4个工具变量",
      global_rssobs = NA_real_, global_p_text = NA_character_, global_p_numeric = NA_real_, global_p_is_bound = NA,
      global_test_significant = NA, raw_beta = NA_real_, raw_se = NA_real_, raw_p = NA_real_,
      corrected_beta = NA_real_, corrected_se = NA_real_, corrected_p = NA_real_,
      outlier_count = NA_integer_, outlier_rsids = NA_character_, distortion_coefficient = NA_real_, distortion_p_text = NA_character_,
      distortion_p_numeric = NA_real_, distortion_p_is_bound = NA, warning_messages = NA_character_, error_message = NA_character_
    )
    next
  }
  if (z$run_status == "ERROR") {
    summary_rows[[i]] <- data.table(
      bac = z$bac, nsnp = z$nsnp, nb_distribution = z$nb_distribution, seed = z$seed,
      run_status = "ERROR", run_status_zh = "MR-PRESSO运行失败",
      global_rssobs = NA_real_, global_p_text = NA_character_, global_p_numeric = NA_real_, global_p_is_bound = NA,
      global_test_significant = NA, raw_beta = NA_real_, raw_se = NA_real_, raw_p = NA_real_,
      corrected_beta = NA_real_, corrected_se = NA_real_, corrected_p = NA_real_,
      outlier_count = NA_integer_, outlier_rsids = NA_character_, distortion_coefficient = NA_real_, distortion_p_text = NA_character_,
      distortion_p_numeric = NA_real_, distortion_p_is_bound = NA, warning_messages = paste(z$warnings, collapse = " | "), error_message = z$error
    )
    next
  }

  result <- z$result
  main <- result[["Main MR results"]]
  presso <- result[["MR-PRESSO results"]]
  global <- presso[["Global Test"]]
  outlier_test <- presso[["Outlier Test"]]
  distortion <- presso[["Distortion Test"]]
  global_p <- parse_p(global$Pvalue)
  raw <- main_row(main, "Raw")
  corrected <- main_row(main, "Outlier-corrected")

  outlier_rsids <- character()
  if (is.data.frame(outlier_test) && nrow(outlier_test) > 0L) {
    labels <- rownames(outlier_test)
    numeric_labels <- suppressWarnings(as.integer(labels))
    labels[!is.na(numeric_labels) & numeric_labels >= 1L & numeric_labels <= length(z$rsids)] <- z$rsids[numeric_labels[!is.na(numeric_labels) & numeric_labels >= 1L & numeric_labels <= length(z$rsids)]]
    for (j in seq_len(nrow(outlier_test))) {
      p <- parse_p(outlier_test[["Pvalue"]][[j]])
      if (!is.finite(p$numeric) || p$numeric > significance_threshold) next
      outlier_rsids <- unique(c(outlier_rsids, labels[[j]]))
      outlier_rows[[outlier_position]] <- data.table(
        bac = z$bac, outlier_rsid = labels[[j]], rssobs = as.numeric(outlier_test[["RSSobs"]][[j]]),
        p_text = p$text, p_numeric = p$numeric, p_is_bound = p$is_bound
      )
      outlier_position <- outlier_position + 1L
    }
  }
  if (!is.null(distortion) && !is.null(distortion[["Outliers Indices"]])) {
    indices <- suppressWarnings(as.integer(distortion[["Outliers Indices"]]))
    indices <- indices[is.finite(indices) & indices >= 1L & indices <= length(z$rsids)]
    outlier_rsids <- unique(c(outlier_rsids, z$rsids[indices]))
  }
  distortion_p <- parse_p(if (is.null(distortion)) NULL else distortion$Pvalue)
  global_sig <- is.finite(global_p$numeric) && global_p$numeric < significance_threshold
  final_status <- if (!global_sig) {
    "OK_NO_GLOBAL_PLEIOTROPY"
  } else if (length(outlier_rsids) > 0L) {
    "OK_GLOBAL_WITH_OUTLIERS"
  } else {
    "OK_GLOBAL_NO_OUTLIERS_IDENTIFIED"
  }
  final_status_zh <- switch(
    final_status,
    OK_NO_GLOBAL_PLEIOTROPY = "完成：MR-PRESSO全局检验未见显著水平多效性",
    OK_GLOBAL_WITH_OUTLIERS = "完成：全局检验显著并识别到异常工具变量",
    OK_GLOBAL_NO_OUTLIERS_IDENTIFIED = "完成：全局检验显著，但未识别到可校正异常工具变量"
  )

  summary_rows[[i]] <- data.table(
    bac = z$bac, nsnp = z$nsnp, nb_distribution = z$nb_distribution, seed = z$seed,
    run_status = final_status, run_status_zh = final_status_zh,
    global_rssobs = as.numeric(global$RSSobs), global_p_text = global_p$text, global_p_numeric = global_p$numeric,
    global_p_is_bound = global_p$is_bound, global_test_significant = global_sig,
    raw_beta = if (is.null(raw)) NA_real_ else as.numeric(raw[["Causal Estimate"]]),
    raw_se = if (is.null(raw)) NA_real_ else as.numeric(raw[["Sd"]]),
    raw_p = if (is.null(raw)) NA_real_ else as.numeric(raw[["P-value"]]),
    corrected_beta = if (is.null(corrected)) NA_real_ else as.numeric(corrected[["Causal Estimate"]]),
    corrected_se = if (is.null(corrected)) NA_real_ else as.numeric(corrected[["Sd"]]),
    corrected_p = if (is.null(corrected)) NA_real_ else as.numeric(corrected[["P-value"]]),
    outlier_count = length(outlier_rsids), outlier_rsids = if (length(outlier_rsids) == 0L) NA_character_ else paste(outlier_rsids, collapse = ";"),
    distortion_coefficient = if (is.null(distortion)) NA_real_ else as.numeric(distortion[["Distortion Coefficient"]]),
    distortion_p_text = distortion_p$text, distortion_p_numeric = distortion_p$numeric, distortion_p_is_bound = distortion_p$is_bound,
    warning_messages = if (length(z$warnings) == 0L) NA_character_ else paste(z$warnings, collapse = " | "), error_message = NA_character_
  )
}

summary_table <- rbindlist(summary_rows, use.names = TRUE, fill = TRUE)
setorder(summary_table, bac)
outlier_table <- if (length(outlier_rows) == 0L) {
  data.table(bac = character(), outlier_rsid = character(), rssobs = numeric(), p_text = character(), p_numeric = numeric(), p_is_bound = logical())
} else {
  rbindlist(outlier_rows, use.names = TRUE, fill = TRUE)
}
data.table::fwrite(summary_table, file.path(output_dir, "stage2_mr_presso_summary.csv"), na = "")
data.table::fwrite(outlier_table, file.path(output_dir, "stage2_mr_presso_outliers.csv"), na = "")

metrics <- list(
  status = if (any(summary_table$run_status == "ERROR")) "MR_PRESSO_COMPLETED_WITH_ERRORS" else "PASS_MR_PRESSO_EXECUTION",
  status_zh = if (any(summary_table$run_status == "ERROR")) "MR-PRESSO完成但存在运行错误，需审计" else "通过：所有适用taxa均完成MR-PRESSO",
  generated_at_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
  specification = list(
    package = "MRPRESSO", package_version = as.character(utils::packageVersion("MRPRESSO")),
    nb_distribution = nb_distribution, significance_threshold = significance_threshold,
    seed_base = seed_base, requested_workers = n_workers
  ),
  results = list(
    taxa = nrow(summary_table),
    applicable_taxa = sum(summary_table$run_status != "NOT_APPLICABLE_NOT_ENOUGH_IVS"),
    not_applicable_taxa = sum(summary_table$run_status == "NOT_APPLICABLE_NOT_ENOUGH_IVS"),
    error_taxa = sum(summary_table$run_status == "ERROR"),
    global_test_significant_taxa = sum(summary_table$global_test_significant %in% TRUE),
    taxa_with_outliers = sum(summary_table$outlier_count > 0, na.rm = TRUE),
    outlier_rows = nrow(outlier_table)
  ),
  software = list(R = R.version.string, MRPRESSO = as.character(utils::packageVersion("MRPRESSO")))
)
jsonlite::write_json(metrics, file.path(qa_dir, "stage2_mr_presso_metrics.json"), pretty = TRUE, auto_unbox = TRUE, na = "null")

cat(sprintf(
  "%s taxa=%d applicable=%d not_applicable=%d errors=%d global_significant=%d taxa_with_outliers=%d\n",
  metrics$status, metrics$results$taxa, metrics$results$applicable_taxa,
  metrics$results$not_applicable_taxa, metrics$results$error_taxa,
  metrics$results$global_test_significant_taxa, metrics$results$taxa_with_outliers
))
