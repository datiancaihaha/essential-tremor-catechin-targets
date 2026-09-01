args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5) {
  stop(paste(
    "Usage: Rscript composition_analysis.R",
    "<cell_counts.csv> <propeller_out> <matched_out> <speckle_library> <analysis_library>"
  ))
}

cell_counts_path <- normalizePath(args[[1]], winslash = "/", mustWork = TRUE)
propeller_out <- args[[2]]
matched_out <- args[[3]]
speckle_library <- normalizePath(args[[4]], winslash = "/", mustWork = TRUE)
analysis_library <- normalizePath(args[[5]], winslash = "/", mustWork = TRUE)
.libPaths(c(analysis_library, speckle_library, .libPaths()))

suppressPackageStartupMessages(library(limma))
suppressPackageStartupMessages(library(speckle))
suppressPackageStartupMessages(library(MatchIt))
options(warn = 1)
dir.create(propeller_out, recursive = TRUE, showWarnings = FALSE)
dir.create(matched_out, recursive = TRUE, showWarnings = FALSE)

cell_counts <- read.csv(
  cell_counts_path,
  check.names = FALSE,
  colClasses = c(donor_id = "character", sample_id = "character")
)
required_columns <- c(
  "donor_id", "cell_type", "seq_batch", "cell_count", "disease",
  "sex", "age", "sample_id"
)
if (!all(required_columns %in% names(cell_counts))) {
  stop("The cell-count input is missing one or more required columns")
}
if (any(cell_counts$cell_count < 0) || any(cell_counts$cell_count %% 1 != 0)) {
  stop("Cell counts must be non-negative integers")
}

cell_types <- sort(unique(cell_counts$cell_type))
if (length(cell_types) != 14) {
  stop(paste("Expected 14 author-defined cell types; found", length(cell_types)))
}

sample_metadata <- unique(cell_counts[, c(
  "sample_id", "donor_id", "disease", "sex", "age", "seq_batch"
)])
if (anyDuplicated(sample_metadata$sample_id)) {
  stop("Sample metadata are not unique by donor-by-sequencing-batch sample")
}
sample_metadata <- sample_metadata[order(sample_metadata$donor_id, sample_metadata$seq_batch), ]

make_count_matrix <- function(rows, column_ids, column_name) {
  count_matrix <- matrix(
    0L,
    nrow = length(cell_types),
    ncol = length(column_ids),
    dimnames = list(cell_types, column_ids)
  )
  row_index <- match(rows$cell_type, cell_types)
  column_index <- match(rows[[column_name]], column_ids)
  if (anyNA(row_index) || anyNA(column_index)) stop("Count-matrix indexing failed")
  for (index in seq_len(nrow(rows))) {
    count_matrix[row_index[[index]], column_index[[index]]] <-
      count_matrix[row_index[[index]], column_index[[index]]] + rows$cell_count[[index]]
  }
  count_matrix
}

sample_count_matrix <- make_count_matrix(
  cell_counts,
  sample_metadata$sample_id,
  "sample_id"
)
observed_sample_totals <- aggregate(
  cell_count ~ sample_id,
  data = cell_counts,
  FUN = sum
)
if (!all(
  colSums(sample_count_matrix) ==
    observed_sample_totals$cell_count[match(colnames(sample_count_matrix), observed_sample_totals$sample_id)]
)) {
  stop("Donor-by-batch cell-count totals are not conserved")
}

donor_metadata <- unique(cell_counts[, c("donor_id", "disease", "sex", "age")])
if (anyDuplicated(donor_metadata$donor_id)) {
  stop("Donor metadata are not unique")
}
donor_metadata <- donor_metadata[order(donor_metadata$donor_id), ]
batch_profiles <- aggregate(
  seq_batch ~ donor_id,
  data = unique(cell_counts[, c("donor_id", "seq_batch")]),
  FUN = function(values) paste(sort(unique(values)), collapse = "+")
)
names(batch_profiles)[2] <- "batch_profile"
donor_metadata <- merge(
  donor_metadata,
  batch_profiles,
  by = "donor_id",
  all.x = TRUE,
  sort = FALSE
)
donor_metadata <- donor_metadata[match(sort(donor_metadata$donor_id), donor_metadata$donor_id), ]
donor_count_matrix <- make_count_matrix(
  cell_counts,
  donor_metadata$donor_id,
  "donor_id"
)
if (sum(donor_count_matrix) != sum(sample_count_matrix)) {
  stop("Donor aggregation did not conserve the total number of cells")
}

add_model_variables <- function(metadata) {
  metadata$condition <- factor(
    ifelse(metadata$disease == "essential tremor", "case", "control"),
    levels = c("control", "case")
  )
  metadata$age_c <- metadata$age - mean(metadata$age)
  metadata$sex_factor <- factor(metadata$sex, levels = c("female", "male"))
  metadata
}

fit_proportions <- function(counts, metadata, model_type) {
  metadata <- metadata[match(colnames(counts), metadata[[
    if (model_type == "donor_total") "donor_id" else "sample_id"
  ]]), , drop = FALSE]
  if (anyNA(metadata$donor_id)) stop("Composition metadata alignment failed")
  metadata <- add_model_variables(metadata)

  if (model_type == "donor_total") {
    metadata$batch_profile_factor <- factor(metadata$batch_profile)
    design <- model.matrix(
      ~ condition + age_c + sex_factor + batch_profile_factor,
      data = metadata
    )
    block <- NULL
  } else {
    metadata$seq_batch_factor <- factor(metadata$seq_batch)
    design <- model.matrix(
      ~ condition + age_c + sex_factor + seq_batch_factor,
      data = metadata
    )
    block <- metadata$donor_id
  }
  if (qr(design)$rank != ncol(design)) {
    stop(paste("Composition design matrix is not full rank for", model_type))
  }

  prop_list <- convertDataToList(counts, data.type = "counts", transform = "logit")
  if (is.null(block)) {
    fit <- lmFit(prop_list$TransformedProps, design)
    consensus_correlation <- NA_real_
  } else {
    correlation <- duplicateCorrelation(
      prop_list$TransformedProps,
      design,
      block = block
    )
    if (!is.finite(correlation$consensus)) {
      stop(paste("Non-finite repeated-measure correlation for", model_type))
    }
    consensus_correlation <- correlation$consensus
    fit <- lmFit(
      prop_list$TransformedProps,
      design,
      block = block,
      correlation = consensus_correlation
    )
  }
  fit <- eBayes(fit, robust = TRUE, trend = FALSE)
  result <- topTable(
    fit,
    coef = "conditioncase",
    number = Inf,
    sort.by = "none",
    adjust.method = "none",
    confint = 0.95
  )
  result$cell_type <- rownames(result)
  rownames(result) <- NULL
  result$BH_FDR_supplementary <- p.adjust(result$P.Value, method = "BH")
  result$nominal_p_lt_0_05 <- result$P.Value < 0.05
  result$consensus_correlation <- consensus_correlation
  result$n_composition_samples <- ncol(counts)
  result$n_donors <- length(unique(metadata$donor_id))
  result$n_et_donors <- length(unique(metadata$donor_id[metadata$condition == "case"]))
  result$n_control_donors <- length(unique(metadata$donor_id[metadata$condition == "control"]))
  result$model_type <- model_type

  donor_means <- donor_count_matrix[, unique(metadata$donor_id), drop = FALSE]
  donor_props <- sweep(donor_means, 2, colSums(donor_means), "/")
  donor_conditions <- donor_metadata$disease[
    match(colnames(donor_props), donor_metadata$donor_id)
  ]
  result$mean_proportion_control <- rowMeans(
    donor_props[, donor_conditions == "normal", drop = FALSE]
  )[result$cell_type]
  result$mean_proportion_et <- rowMeans(
    donor_props[, donor_conditions == "essential tremor", drop = FALSE]
  )[result$cell_type]
  result$mean_proportion_difference_et_minus_control <-
    result$mean_proportion_et - result$mean_proportion_control
  result$mean_proportion_ratio_et_over_control <-
    result$mean_proportion_et / result$mean_proportion_control

  preferred_columns <- c(
    "cell_type", "logFC", "CI.L", "CI.R", "t", "P.Value",
    "nominal_p_lt_0_05", "BH_FDR_supplementary",
    "mean_proportion_control", "mean_proportion_et",
    "mean_proportion_difference_et_minus_control",
    "mean_proportion_ratio_et_over_control", "consensus_correlation",
    "n_composition_samples", "n_donors", "n_et_donors",
    "n_control_donors", "model_type"
  )
  list(
    result = result[, preferred_columns],
    design = design,
    proportions = prop_list$Proportions,
    metadata = metadata
  )
}

main_fit <- fit_proportions(
  sample_count_matrix,
  sample_metadata,
  "donor_by_sequencing_batch_blocked"
)
main_result <- main_fit$result[order(main_fit$result$P.Value), ]
write.csv(
  main_result,
  file.path(propeller_out, "propeller_blocked_all_14_cell_types.csv"),
  row.names = FALSE
)
write.csv(
  main_result[main_result$nominal_p_lt_0_05, ],
  file.path(propeller_out, "propeller_blocked_nominal_candidates.csv"),
  row.names = FALSE
)

donor_total_fit <- fit_proportions(
  donor_count_matrix,
  donor_metadata,
  "donor_total"
)
donor_total_result <- donor_total_fit$result[order(donor_total_fit$result$P.Value), ]
write.csv(
  donor_total_result,
  file.path(propeller_out, "propeller_donor_total_all_14_cell_types.csv"),
  row.names = FALSE
)

nominal_cell_types <- main_result$cell_type[main_result$nominal_p_lt_0_05]
et_donors <- donor_metadata$donor_id[donor_metadata$disease == "essential tremor"]
leave_one_out <- list()
if (length(nominal_cell_types) > 0) {
  for (omitted_donor in et_donors) {
    keep_samples <- sample_metadata$donor_id != omitted_donor
    fit <- fit_proportions(
      sample_count_matrix[, keep_samples, drop = FALSE],
      sample_metadata[keep_samples, , drop = FALSE],
      "donor_by_sequencing_batch_blocked"
    )
    result <- fit$result[
      match(nominal_cell_types, fit$result$cell_type),
      c("cell_type", "logFC", "P.Value", "nominal_p_lt_0_05"),
      drop = FALSE
    ]
    result$omitted_et_donor <- omitted_donor
    result$main_logFC <- main_result$logFC[
      match(result$cell_type, main_result$cell_type)
    ]
    result$same_direction_as_main <- sign(result$logFC) == sign(result$main_logFC)
    leave_one_out[[omitted_donor]] <- result
  }
}
if (length(leave_one_out) > 0) {
  leave_one_out <- do.call(rbind, leave_one_out)
  rownames(leave_one_out) <- NULL
  write.csv(
    leave_one_out,
    file.path(propeller_out, "leave_one_et_donor_out_results.csv"),
    row.names = FALSE
  )
  leave_one_out_summary <- do.call(rbind, lapply(
    split(leave_one_out, leave_one_out$cell_type),
    function(values) {
      data.frame(
        cell_type = values$cell_type[[1]],
        omitted_et_donors = nrow(values),
        same_direction_fraction = mean(values$same_direction_as_main),
        nominal_p_lt_0_05_fraction = mean(values$nominal_p_lt_0_05),
        maximum_nominal_p = max(values$P.Value),
        minimum_logFC = min(values$logFC),
        maximum_logFC = max(values$logFC)
      )
    }
  ))
  rownames(leave_one_out_summary) <- NULL
  write.csv(
    leave_one_out_summary,
    file.path(propeller_out, "leave_one_et_donor_out_summary.csv"),
    row.names = FALSE
  )
}

donor_for_matching <- add_model_variables(donor_metadata)
match_object <- matchit(
  condition ~ age,
  data = donor_for_matching,
  method = "nearest",
  distance = "mahalanobis",
  exact = ~ sex_factor,
  ratio = 1,
  replace = FALSE,
  estimand = "ATT"
)
matched_data <- match.data(match_object, data = donor_for_matching)
matched_donors <- matched_data$donor_id
if (
  sum(matched_data$condition == "case") != 16 ||
  sum(matched_data$condition == "control") != 16 ||
  anyDuplicated(matched_donors)
) {
  stop("Age- and sex-matched sensitivity did not retain 16 unique donors per group")
}
write.csv(
  matched_data[, c(
    "donor_id", "disease", "condition", "sex", "age", "batch_profile",
    "weights", "subclass"
  )],
  file.path(matched_out, "matched_donor_metadata.csv"),
  row.names = FALSE
)

standardized_difference <- function(values, groups) {
  case_values <- values[groups == "case"]
  control_values <- values[groups == "control"]
  pooled_sd <- sqrt((var(case_values) + var(control_values)) / 2)
  if (pooled_sd == 0) return(0)
  (mean(case_values) - mean(control_values)) / pooled_sd
}
balance_table <- rbind(
  data.frame(
    covariate = "age",
    standardized_difference_before = standardized_difference(
      donor_for_matching$age,
      donor_for_matching$condition
    ),
    standardized_difference_after = standardized_difference(
      matched_data$age,
      matched_data$condition
    )
  ),
  data.frame(
    covariate = "male",
    standardized_difference_before = standardized_difference(
      as.integer(donor_for_matching$sex == "male"),
      donor_for_matching$condition
    ),
    standardized_difference_after = standardized_difference(
      as.integer(matched_data$sex == "male"),
      matched_data$condition
    )
  )
)
write.csv(
  balance_table,
  file.path(matched_out, "matching_balance.csv"),
  row.names = FALSE
)

matched_sample_metadata <- sample_metadata[
  sample_metadata$donor_id %in% matched_donors,
  ,
  drop = FALSE
]
matched_sample_count_matrix <- sample_count_matrix[
  , matched_sample_metadata$sample_id, drop = FALSE
]
matched_fit <- fit_proportions(
  matched_sample_count_matrix,
  matched_sample_metadata,
  "donor_by_sequencing_batch_blocked"
)
matched_result <- matched_fit$result[order(matched_fit$result$P.Value), ]
write.csv(
  matched_result,
  file.path(matched_out, "matched_propeller_all_14_cell_types.csv"),
  row.names = FALSE
)
write.csv(
  matched_result[matched_result$nominal_p_lt_0_05, ],
  file.path(matched_out, "matched_propeller_nominal_candidates.csv"),
  row.names = FALSE
)

composition_comparison <- merge(
  main_result[, c(
    "cell_type", "logFC", "P.Value", "nominal_p_lt_0_05"
  )],
  donor_total_result[, c(
    "cell_type", "logFC", "P.Value", "nominal_p_lt_0_05"
  )],
  by = "cell_type",
  suffixes = c("_blocked_main", "_donor_total")
)
composition_comparison <- merge(
  composition_comparison,
  matched_result[, c(
    "cell_type", "logFC", "P.Value", "nominal_p_lt_0_05"
  )],
  by = "cell_type"
)
names(composition_comparison)[
  names(composition_comparison) %in% c("logFC", "P.Value", "nominal_p_lt_0_05")
] <- paste0(
  names(composition_comparison)[
    names(composition_comparison) %in% c("logFC", "P.Value", "nominal_p_lt_0_05")
  ],
  "_matched"
)
composition_comparison$same_direction_all_models <- with(
  composition_comparison,
  sign(logFC_blocked_main) == sign(logFC_donor_total) &
    sign(logFC_blocked_main) == sign(logFC_matched)
)
write.csv(
  composition_comparison,
  file.path(propeller_out, "composition_model_comparison.csv"),
  row.names = FALSE
)

sccoda_input <- donor_metadata[, c(
  "donor_id", "disease", "sex", "age", "batch_profile"
)]
sccoda_input$condition <- ifelse(
  sccoda_input$disease == "essential tremor", "case", "control"
)
sccoda_input$age_c <- sccoda_input$age - mean(sccoda_input$age)
for (cell_type in cell_types) {
  sccoda_input[[cell_type]] <- donor_count_matrix[
    cell_type,
    match(sccoda_input$donor_id, colnames(donor_count_matrix))
  ]
}
write.csv(
  sccoda_input,
  file.path(propeller_out, "sccoda_donor_counts_input.csv"),
  row.names = FALSE
)
write.csv(
  sccoda_input[sccoda_input$donor_id %in% matched_donors, ],
  file.path(matched_out, "sccoda_matched_donor_counts_input.csv"),
  row.names = FALSE
)

design_audit <- rbind(
  data.frame(
    model = "blocked_main",
    rows = nrow(main_fit$design),
    columns = ncol(main_fit$design),
    rank = qr(main_fit$design)$rank
  ),
  data.frame(
    model = "donor_total",
    rows = nrow(donor_total_fit$design),
    columns = ncol(donor_total_fit$design),
    rank = qr(donor_total_fit$design)$rank
  ),
  data.frame(
    model = "matched_blocked",
    rows = nrow(matched_fit$design),
    columns = ncol(matched_fit$design),
    rank = qr(matched_fit$design)$rank
  )
)
write.csv(
  design_audit,
  file.path(propeller_out, "composition_design_audit.csv"),
  row.names = FALSE
)

summary_lines <- c(
  "# SCP3177细胞组成分析",
  "",
  "## 分析范围",
  "",
  paste0(
    "纳入", nrow(donor_metadata), "名供体和", length(cell_types),
    "类作者定义细胞类型；ET供体", sum(donor_metadata$disease == "essential tremor"),
    "名，对照供体", sum(donor_metadata$disease == "normal"), "名。"
  ),
  "主模型对供体×测序批次的细胞比例作propeller logit变换，以ET状态、年龄、性别和测序批次为协变量，并以供体作为重复测量块。",
  "另行完成供体单一聚合模型，以及16名ET供体与16名年龄和性别匹配对照的敏感性分析。",
  "候选细胞类型由nominal P < 0.05定义；BH校正值仅作为补充证据强度报告。",
  "",
  "## 完成度",
  "",
  paste0("主模型14类细胞均获得完整结果；", sum(main_result$nominal_p_lt_0_05), "类达到nominal P < 0.05。"),
  paste0("匹配对照模型中，", sum(matched_result$nominal_p_lt_0_05), "类达到nominal P < 0.05。"),
  "",
  "## 方法依据",
  "",
  "- Phipson B, et al. Bioinformatics. 2022;38:4720-4726. doi:10.1093/bioinformatics/btac582.",
  "- Castonguay C, et al. Nature Genetics. 2026. doi:10.1038/s41588-026-02544-8.",
  "- Rosenbaum PR. Annual Review of Statistics and Its Application. 2020;7:143-176. doi:10.1146/annurev-statistics-031219-041058."
)
writeLines(
  summary_lines,
  file.path(propeller_out, "CELL_COMPOSITION_RESULTS_zh.md"),
  useBytes = TRUE
)
capture.output(sessionInfo(), file = file.path(propeller_out, "session_info.txt"))
