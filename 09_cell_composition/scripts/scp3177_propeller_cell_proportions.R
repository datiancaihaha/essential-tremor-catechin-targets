args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
  stop("Usage: Rscript scp3177_propeller_cell_proportions.R <cell_counts.csv> <out_dir> <r_library>")
}

cell_counts_path <- normalizePath(args[[1]], winslash = "/", mustWork = TRUE)
out_dir <- args[[2]]
r_library <- normalizePath(args[[3]], winslash = "/", mustWork = TRUE)
.libPaths(c(r_library, .libPaths()))

suppressPackageStartupMessages(library(limma))
suppressPackageStartupMessages(library(speckle))
options(warn = 1)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

counts_long <- read.csv(cell_counts_path, check.names = FALSE)
counts_long$donor_id <- as.character(counts_long$donor_id)
donor_metadata <- unique(counts_long[, c("donor_id", "disease", "sex", "age", "total_cells")])
donor_metadata <- donor_metadata[order(donor_metadata$donor_id), ]

cell_types <- sort(unique(counts_long$cell_type))
donors <- donor_metadata$donor_id
count_matrix <- matrix(
  0,
  nrow = length(cell_types),
  ncol = length(donors),
  dimnames = list(cell_types, donors)
)
count_matrix[cbind(
  match(counts_long$cell_type, cell_types),
  match(counts_long$donor_id, donors)
)] <- counts_long$cell_count

if (!all(colSums(count_matrix) == donor_metadata$total_cells)) {
  stop("Cell-type counts do not sum to donor total_cells")
}

run_propeller <- function(counts, metadata) {
  metadata <- metadata[match(colnames(counts), metadata$donor_id), , drop = FALSE]
  if (anyNA(metadata$donor_id)) stop("Donor metadata alignment failed")
  metadata$group <- factor(
    metadata$disease,
    levels = c("normal", "essential tremor"),
    labels = c("Normal", "ET")
  )
  metadata$age_c <- metadata$age - mean(metadata$age)
  metadata$male <- as.integer(metadata$sex == "male")
  design <- model.matrix(~ 0 + group + age_c + male, data = metadata)
  if (qr(design)$rank != ncol(design)) stop("Cell-proportion design matrix is not full rank")
  contrast <- makeContrasts(ETvsNormal = groupET - groupNormal, levels = design)
  prop_list <- convertDataToList(counts, data.type = "counts", transform = "logit")
  result <- propeller.ttest(
    prop.list = prop_list,
    design = design,
    contrasts = contrast[, "ETvsNormal"],
    robust = TRUE,
    trend = FALSE,
    sort = FALSE
  )
  list(result = result, proportions = prop_list$Proportions)
}

main <- run_propeller(count_matrix, donor_metadata)
results <- main$result
results$cell_type <- rownames(results)
rownames(results) <- NULL
results$mean_proportion_normal <- rowMeans(
  main$proportions[, donor_metadata$disease == "normal", drop = FALSE]
)
results$mean_proportion_et <- rowMeans(
  main$proportions[, donor_metadata$disease == "essential tremor", drop = FALSE]
)
results$mean_proportion_difference_et_minus_normal <-
  results$mean_proportion_et - results$mean_proportion_normal
results$nominal_p_lt_0_05 <- results$P.Value < 0.05
names(results)[names(results) == "FDR"] <- "BH_FDR_supplementary"
results$n_donors <- nrow(donor_metadata)
results$n_et <- sum(donor_metadata$disease == "essential tremor")
results$n_normal <- sum(donor_metadata$disease == "normal")
results$analysis_type <- "donor-level cell-type proportion analysis"

preferred_columns <- c(
  "cell_type", "mean_proportion_normal", "mean_proportion_et",
  "mean_proportion_difference_et_minus_normal", "PropRatio", "Tstatistic",
  "P.Value", "nominal_p_lt_0_05", "BH_FDR_supplementary",
  "n_donors", "n_et", "n_normal", "analysis_type"
)
results <- results[, preferred_columns]
results <- results[order(results$P.Value), ]
write.csv(results, file.path(out_dir, "propeller_all_cell_types.csv"), row.names = FALSE)
write.csv(
  results[results$nominal_p_lt_0_05, ],
  file.path(out_dir, "propeller_nominal_candidates.csv"),
  row.names = FALSE
)

proportions_long <- counts_long[, c(
  "donor_id", "cell_type", "cell_count", "disease", "sex", "age", "total_cells"
)]
proportions_long$cell_fraction <- proportions_long$cell_count / proportions_long$total_cells
write.csv(
  proportions_long,
  file.path(out_dir, "donor_cell_type_proportions.csv"),
  row.names = FALSE
)

nominal_cell_types <- results$cell_type[results$nominal_p_lt_0_05]
case_donors <- donor_metadata$donor_id[donor_metadata$disease == "essential tremor"]
loo_results <- list()
if (length(nominal_cell_types) > 0) {
  for (omitted_donor in case_donors) {
    keep <- colnames(count_matrix) != omitted_donor
    metadata_loo <- donor_metadata[donor_metadata$donor_id != omitted_donor, , drop = FALSE]
    loo <- run_propeller(
      count_matrix[, keep, drop = FALSE],
      metadata_loo
    )
    loo_result <- loo$result[nominal_cell_types, , drop = FALSE]
    loo_result$cell_type <- rownames(loo_result)
    loo_result$omitted_et_donor <- omitted_donor
    loo_mean_normal <- rowMeans(
      loo$proportions[, metadata_loo$disease == "normal", drop = FALSE]
    )
    loo_mean_et <- rowMeans(
      loo$proportions[, metadata_loo$disease == "essential tremor", drop = FALSE]
    )
    loo_result$mean_proportion_normal <- loo_mean_normal[loo_result$cell_type]
    loo_result$mean_proportion_et <- loo_mean_et[loo_result$cell_type]
    loo_result$mean_proportion_difference_et_minus_normal <-
      loo_result$mean_proportion_et - loo_result$mean_proportion_normal
    loo_result$nominal_p_lt_0_05 <- loo_result$P.Value < 0.05
    loo_results[[omitted_donor]] <- loo_result[, c(
      "cell_type", "omitted_et_donor", "mean_proportion_difference_et_minus_normal",
      "PropRatio", "Tstatistic", "P.Value", "nominal_p_lt_0_05"
    )]
  }
}

if (length(loo_results) > 0) {
  loo <- do.call(rbind, loo_results)
  rownames(loo) <- NULL
  main_direction <- results[, c("cell_type", "mean_proportion_difference_et_minus_normal")]
  names(main_direction)[2] <- "main_mean_proportion_difference"
  loo <- merge(loo, main_direction, by = "cell_type", all.x = TRUE)
  loo$same_direction_as_main <- sign(loo$mean_proportion_difference_et_minus_normal) ==
    sign(loo$main_mean_proportion_difference)
  write.csv(
    loo,
    file.path(out_dir, "leave_one_et_donor_out_results.csv"),
    row.names = FALSE
  )
  loo_summary <- do.call(rbind, lapply(split(loo, loo$cell_type), function(values) {
    data.frame(
      cell_type = values$cell_type[[1]],
      omitted_et_donors = nrow(values),
      same_direction_fraction = mean(values$same_direction_as_main),
      nominal_p_lt_0_05_fraction = mean(values$nominal_p_lt_0_05),
      maximum_nominal_p = max(values$P.Value),
      minimum_mean_proportion_difference = min(values$mean_proportion_difference_et_minus_normal),
      maximum_mean_proportion_difference = max(values$mean_proportion_difference_et_minus_normal)
    )
  }))
  rownames(loo_summary) <- NULL
  write.csv(
    loo_summary,
    file.path(out_dir, "leave_one_et_donor_out_summary.csv"),
    row.names = FALSE
  )
}

summary_lines <- c(
  "# SCP3177供体层细胞组成分析",
  "",
  "## 方法",
  "",
  paste0(
    "以", nrow(donor_metadata), "名供体作为生物学重复（ET ",
    sum(donor_metadata$disease == "essential tremor"), "名；对照 ",
    sum(donor_metadata$disease == "normal"), "名）。"
  ),
  paste0(
    "使用speckle ", as.character(packageVersion("speckle")),
    "的propeller方法，对13类细胞的供体内比例作logit变换并采用0.5伪计数，",
    "线性模型校正中心化年龄和性别，稳健经验贝叶斯估计组间差异。"
  ),
  "候选纳入、保留和展示由nominal P < 0.05决定；BH校正值仅作为补充证据强度记录。",
  "对主分析nominal P < 0.05的细胞类型逐一剔除每名ET供体，评价方向和名义显著性的稳定性。",
  "",
  "## 结果规模",
  "",
  paste0("13类细胞中，", sum(results$nominal_p_lt_0_05), "类达到nominal P < 0.05。"),
  paste0(
    "其中", sum(results$BH_FDR_supplementary < 0.05),
    "类的补充BH校正值小于0.05。"
  ),
  "",
  "## 方法学文献",
  "",
  "- Phipson et al. Bioinformatics. 2022. doi:10.1093/bioinformatics/btac582.",
  "- Castonguay et al. Nature Genetics. 2026. doi:10.1038/s41588-026-02544-8."
)
writeLines(
  summary_lines,
  file.path(out_dir, "CELL_PROPORTION_ANALYSIS_SUMMARY_zh.md"),
  useBytes = TRUE
)
capture.output(sessionInfo(), file = file.path(out_dir, "session_info.txt"))
