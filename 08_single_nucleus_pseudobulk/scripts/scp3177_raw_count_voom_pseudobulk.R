args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
  stop("Usage: Rscript scp3177_raw_count_voom_pseudobulk.R <pseudobulk_dir> <target_matrix.csv> <out_dir>")
}

suppressPackageStartupMessages(library(edgeR))
suppressPackageStartupMessages(library(limma))
options(warn = 1)

pseudobulk_dir <- normalizePath(args[[1]], winslash = "/", mustWork = TRUE)
target_path <- normalizePath(args[[2]], winslash = "/", mustWork = TRUE)
out_dir <- args[[3]]
full_result_dir <- file.path(out_dir, "full_transcriptome_results")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(full_result_dir, recursive = TRUE, showWarnings = FALSE)

genes <- read.csv(file.path(pseudobulk_dir, "genes.csv"), check.names = FALSE)
manifest <- read.csv(
  file.path(pseudobulk_dir, "aggregation_manifest.csv"),
  check.names = FALSE
)
targets <- read.csv(
  target_path,
  check.names = FALSE,
  fileEncoding = "UTF-8-BOM"
)
prediction_columns <- grep("_predicted$", names(targets), value = TRUE)
targets[prediction_columns] <- lapply(
  targets[prediction_columns],
  function(values) tolower(as.character(values)) == "true"
)
targets$all_form_union <- rowSums(targets[prediction_columns] == TRUE) > 0
targets$all_form_core <- rowSums(targets[prediction_columns] == TRUE) == length(prediction_columns)
circulating_columns <- prediction_columns[
  grepl("SULFATE|GLUCURONIDE", prediction_columns)
]
targets$circulating_conjugate_union <-
  rowSums(targets[circulating_columns] == TRUE) > 0
targets$circulating_conjugate_core <-
  rowSums(targets[circulating_columns] == TRUE) == length(circulating_columns)
prespecified_genes <- c("CA1", "CA2", "CA3", "CA13", "ADCY10", "CA7", "TET3", "SIRT1", "IGF1R")
targets$prespecified_single_gene <- targets$gene_symbol %in% prespecified_genes
target_flags <- targets[, c(
  "gene_symbol", "all_form_union", "all_form_core",
  "circulating_conjugate_union", "circulating_conjugate_core",
  "prespecified_single_gene"
)]
if (anyDuplicated(target_flags$gene_symbol)) stop("Duplicate target gene symbols")
if (!all(target_flags$gene_symbol %in% genes$gene_symbol)) {
  stop("One or more target genes are absent from the pseudobulk gene table")
}

all_sample_metadata <- read.csv(
  file.path(pseudobulk_dir, "all_group_metadata.csv"),
  check.names = FALSE,
  colClasses = c(donor_id = "character", sample_id = "character")
)
seq_batch_levels <- sort(unique(all_sample_metadata$seq_batch))

read_count_matrix <- function(slug, n_genes, n_samples) {
  path <- file.path(
    pseudobulk_dir,
    "count_matrices",
    paste0(slug, ".gene_by_sample.int32.bin")
  )
  connection <- file(path, open = "rb")
  on.exit(close(connection))
  values <- readBin(
    connection,
    what = integer(),
    n = n_genes * n_samples,
    size = 4,
    signed = TRUE,
    endian = "little"
  )
  if (length(values) != n_genes * n_samples) {
    stop(paste("Unexpected binary matrix length for", slug))
  }
  matrix(values, nrow = n_genes, ncol = n_samples)
}

prepare_design <- function(metadata) {
  metadata$condition <- factor(
    ifelse(metadata$disease == "essential tremor", "case", "control"),
    levels = c("control", "case")
  )
  metadata$age_c <- metadata$age - mean(metadata$age)
  metadata$sex_factor <- factor(metadata$sex, levels = c("female", "male"))
  metadata$seq_batch_factor <- factor(
    metadata$seq_batch,
    levels = seq_batch_levels
  )
  design <- model.matrix(
    ~ condition + age_c + sex_factor + seq_batch_factor,
    data = metadata
  )
  if (qr(design)$rank != ncol(design)) {
    stop("Pseudobulk design matrix is not full rank")
  }
  list(metadata = metadata, design = design)
}

fit_voom <- function(counts, metadata) {
  prepared <- prepare_design(metadata)
  metadata <- prepared$metadata
  design <- prepared$design
  dge <- DGEList(counts = counts)
  keep <- filterByExpr(dge, group = metadata$condition)
  if (sum(keep) == 0) stop("No genes passed filterByExpr")
  dge <- dge[keep, , keep.lib.sizes = FALSE]
  dge <- calcNormFactors(dge)
  voom_initial <- voom(dge, design, plot = FALSE)
  correlation_initial <- duplicateCorrelation(
    voom_initial,
    design,
    block = metadata$donor_id
  )
  if (!is.finite(correlation_initial$consensus)) {
    stop("Initial duplicateCorrelation estimate is not finite")
  }
  voom_final <- voom(
    dge,
    design,
    plot = FALSE,
    block = metadata$donor_id,
    correlation = correlation_initial$consensus
  )
  correlation_final <- duplicateCorrelation(
    voom_final,
    design,
    block = metadata$donor_id
  )
  if (!is.finite(correlation_final$consensus)) {
    stop("Final duplicateCorrelation estimate is not finite")
  }
  fit <- lmFit(
    voom_final,
    design,
    block = metadata$donor_id,
    correlation = correlation_final$consensus
  )
  fit <- eBayes(fit)
  table <- topTable(
    fit,
    coef = "conditioncase",
    number = Inf,
    sort.by = "none",
    adjust.method = "BH",
    confint = 0.95
  )
  list(
    table = table,
    keep = keep,
    normalization_factors = dge$samples$norm.factors,
    initial_correlation = correlation_initial$consensus,
    final_correlation = correlation_final$consensus,
    retained_genes = sum(keep)
  )
}

main_results <- list()
filter_audit <- list()
model_summary <- list()

for (manifest_index in seq_len(nrow(manifest))) {
  cell_type <- manifest$cell_type[[manifest_index]]
  slug <- manifest$slug[[manifest_index]]
  metadata <- read.csv(
    file.path(pseudobulk_dir, "sample_metadata", paste0(slug, ".csv")),
    check.names = FALSE,
    colClasses = c(donor_id = "character", sample_id = "character")
  )
  counts <- read_count_matrix(slug, nrow(genes), nrow(metadata))
  rownames(counts) <- genes$gene_symbol
  colnames(counts) <- metadata$sample_id
  if (!all(colSums(counts) == metadata$library_size)) {
    stop(paste("Library-size mismatch for", cell_type))
  }

  fitted <- fit_voom(counts, metadata)
  full_table <- fitted$table
  full_table$gene_symbol <- rownames(full_table)
  full_table$gene_id <- genes$gene_id[match(full_table$gene_symbol, genes$gene_symbol)]
  full_table$cell_type <- cell_type
  names(full_table)[names(full_table) == "adj.P.Val"] <-
    "BH_FDR_full_transcriptome_within_cell_type_supplementary"
  full_result_connection <- gzfile(
    file.path(full_result_dir, paste0(slug, ".full_transcriptome.csv.gz")),
    open = "wt"
  )
  write.csv(
    full_table,
    full_result_connection,
    row.names = FALSE
  )
  close(full_result_connection)

  target_table <- full_table[
    full_table$gene_symbol %in% target_flags$gene_symbol,
    ,
    drop = FALSE
  ]
  target_table <- merge(
    target_table,
    target_flags,
    by = "gene_symbol",
    all.x = TRUE,
    sort = FALSE
  )
  target_table$n_pseudobulk_samples <- nrow(metadata)
  target_table$n_donors <- length(unique(metadata$donor_id))
  target_table$n_et_donors <- length(unique(
    metadata$donor_id[metadata$disease == "essential tremor"]
  ))
  target_table$n_normal_donors <- length(unique(
    metadata$donor_id[metadata$disease == "normal"]
  ))
  target_table$initial_consensus_correlation <- fitted$initial_correlation
  target_table$final_consensus_correlation <- fitted$final_correlation
  target_table$nominal_p_lt_0_05 <- target_table$P.Value < 0.05
  target_table$analysis_type <- "raw-count donor-by-sequencing-batch pseudobulk voom analysis"
  main_results[[cell_type]] <- target_table

  target_indices <- match(target_flags$gene_symbol, genes$gene_symbol)
  filter_audit[[cell_type]] <- data.frame(
    cell_type = cell_type,
    gene_symbol = target_flags$gene_symbol,
    total_raw_count = rowSums(counts[target_indices, , drop = FALSE]),
    samples_with_positive_count = rowSums(counts[target_indices, , drop = FALSE] > 0),
    passed_filterByExpr = fitted$keep[target_indices]
  )
  model_summary[[cell_type]] <- data.frame(
    cell_type = cell_type,
    pseudobulk_samples = nrow(metadata),
    donors = length(unique(metadata$donor_id)),
    et_donors = length(unique(metadata$donor_id[metadata$disease == "essential tremor"])),
    normal_donors = length(unique(metadata$donor_id[metadata$disease == "normal"])),
    genes_before_filter = nrow(counts),
    genes_after_filterByExpr = fitted$retained_genes,
    minimum_normalization_factor = min(fitted$normalization_factors),
    maximum_normalization_factor = max(fitted$normalization_factors),
    initial_consensus_correlation = fitted$initial_correlation,
    final_consensus_correlation = fitted$final_correlation
  )
  message(
    "MAIN_MODEL_COMPLETE cell_type=", cell_type,
    " retained_genes=", fitted$retained_genes,
    " target_tests=", nrow(target_table)
  )
  rm(counts, fitted, full_table, target_table)
  invisible(gc(FALSE))
}

all_results <- do.call(rbind, main_results)
rownames(all_results) <- NULL
all_results$BH_FDR_all_target_gene_cell_type_tests_supplementary <- p.adjust(
  all_results$P.Value,
  method = "BH"
)
prespecified_index <- which(all_results$prespecified_single_gene %in% TRUE)
prespecified_q <- rep(NA_real_, nrow(all_results))
prespecified_q[prespecified_index] <- p.adjust(
  all_results$P.Value[prespecified_index],
  method = "BH"
)
all_results$BH_FDR_prespecified_gene_cell_type_tests_supplementary <- prespecified_q

preferred_columns <- c(
  "cell_type", "gene_symbol", "gene_id", "logFC", "CI.L", "CI.R",
  "AveExpr", "t", "P.Value", "nominal_p_lt_0_05",
  "BH_FDR_full_transcriptome_within_cell_type_supplementary",
  "BH_FDR_all_target_gene_cell_type_tests_supplementary",
  "BH_FDR_prespecified_gene_cell_type_tests_supplementary",
  "n_pseudobulk_samples", "n_donors", "n_et_donors", "n_normal_donors",
  "initial_consensus_correlation", "final_consensus_correlation",
  "all_form_union", "all_form_core", "circulating_conjugate_union",
  "circulating_conjugate_core", "prespecified_single_gene", "B", "analysis_type"
)
all_results <- all_results[, preferred_columns]
all_results <- all_results[order(all_results$P.Value), ]
write.csv(
  all_results,
  file.path(out_dir, "all_target_gene_cell_type_results.csv"),
  row.names = FALSE
)
write.csv(
  all_results[all_results$nominal_p_lt_0_05, ],
  file.path(out_dir, "nominal_target_candidates.csv"),
  row.names = FALSE
)
write.csv(
  all_results[all_results$prespecified_single_gene, ],
  file.path(out_dir, "prespecified_gene_results.csv"),
  row.names = FALSE
)
filter_audit <- do.call(rbind, filter_audit)
rownames(filter_audit) <- NULL
write.csv(filter_audit, file.path(out_dir, "target_filter_audit.csv"), row.names = FALSE)
model_summary <- do.call(rbind, model_summary)
rownames(model_summary) <- NULL
write.csv(model_summary, file.path(out_dir, "cell_type_model_summary.csv"), row.names = FALSE)

loo_results <- list()
nominal_by_cell_type <- split(
  all_results$gene_symbol[all_results$nominal_p_lt_0_05],
  all_results$cell_type[all_results$nominal_p_lt_0_05]
)
case_donors <- sort(unique(
  all_sample_metadata$donor_id[
    all_sample_metadata$disease == "essential tremor"
  ]
))

for (manifest_index in seq_len(nrow(manifest))) {
  cell_type <- manifest$cell_type[[manifest_index]]
  if (!cell_type %in% names(nominal_by_cell_type)) next
  nominal_genes <- nominal_by_cell_type[[cell_type]]
  slug <- manifest$slug[[manifest_index]]
  metadata <- read.csv(
    file.path(pseudobulk_dir, "sample_metadata", paste0(slug, ".csv")),
    check.names = FALSE,
    colClasses = c(donor_id = "character", sample_id = "character")
  )
  counts <- read_count_matrix(slug, nrow(genes), nrow(metadata))
  rownames(counts) <- genes$gene_symbol
  colnames(counts) <- metadata$sample_id

  for (donor_index in seq_along(case_donors)) {
    omitted_donor <- case_donors[[donor_index]]
    keep_samples <- metadata$donor_id != omitted_donor
    fitted <- fit_voom(
      counts[, keep_samples, drop = FALSE],
      metadata[keep_samples, , drop = FALSE]
    )
    table <- fitted$table
    table$gene_symbol <- rownames(table)
    matches <- match(nominal_genes, table$gene_symbol)
    main_lookup <- all_results[
      all_results$cell_type == cell_type &
        all_results$gene_symbol %in% nominal_genes,
      c("gene_symbol", "logFC"),
      drop = FALSE
    ]
    main_logfc <- main_lookup$logFC[match(nominal_genes, main_lookup$gene_symbol)]
    loo_results[[paste(cell_type, omitted_donor, sep = "|")]] <- data.frame(
      cell_type = cell_type,
      gene_symbol = nominal_genes,
      omitted_et_donor = omitted_donor,
      passed_filterByExpr = !is.na(matches),
      logFC = ifelse(is.na(matches), NA_real_, table$logFC[matches]),
      P.Value = ifelse(is.na(matches), NA_real_, table$P.Value[matches]),
      main_logFC = main_logfc
    )
    if (donor_index %% 4 == 0 || donor_index == length(case_donors)) {
      message(
        "LOO_PROGRESS cell_type=", cell_type,
        " donors=", donor_index, "/", length(case_donors)
      )
      invisible(gc(FALSE))
    }
  }
  rm(counts)
  invisible(gc(FALSE))
}

loo <- do.call(rbind, loo_results)
rownames(loo) <- NULL
loo$same_direction_as_main <-
  sign(loo$logFC) == sign(loo$main_logFC)
loo$nominal_p_lt_0_05 <- !is.na(loo$P.Value) & loo$P.Value < 0.05
write.csv(
  loo,
  file.path(out_dir, "leave_one_et_donor_out_results.csv"),
  row.names = FALSE
)
loo_summary <- do.call(rbind, lapply(
  split(loo, interaction(loo$cell_type, loo$gene_symbol, drop = TRUE)),
  function(values) {
    data.frame(
      cell_type = values$cell_type[[1]],
      gene_symbol = values$gene_symbol[[1]],
      omitted_et_donors = nrow(values),
      filter_pass_fraction = mean(values$passed_filterByExpr),
      same_direction_fraction = mean(values$same_direction_as_main, na.rm = TRUE),
      nominal_p_lt_0_05_fraction = mean(values$nominal_p_lt_0_05),
      all_iterations_nominal_p_lt_0_05 = all(values$nominal_p_lt_0_05),
      maximum_nominal_p = if (all(is.na(values$P.Value))) NA_real_ else max(values$P.Value, na.rm = TRUE),
      minimum_logFC = if (all(is.na(values$logFC))) NA_real_ else min(values$logFC, na.rm = TRUE),
      maximum_logFC = if (all(is.na(values$logFC))) NA_real_ else max(values$logFC, na.rm = TRUE)
    )
  }
))
rownames(loo_summary) <- NULL
write.csv(
  loo_summary,
  file.path(out_dir, "leave_one_et_donor_out_summary.csv"),
  row.names = FALSE
)

summary_lines <- c(
  "# SCP3177 raw-count pseudobulk分析",
  "",
  "## 方法",
  "",
  paste0(
    "以供体×seq_batch×细胞类型聚合raw counts；每个pseudobulk样本至少包含10个细胞。",
    "共分析", nrow(manifest), "类细胞。"
  ),
  "使用edgeR::filterByExpr和TMM标准化，并按公开作者代码运行两轮voom/duplicateCorrelation。",
  "模型包含ET状态、中心化年龄、性别和seq_batch，供体作为重复测量块。",
  "候选纳入、保留和展示由nominal P < 0.05决定；BH校正值仅作为补充证据强度记录。",
  "",
  "## 结果规模",
  "",
  paste0(
    "完成", nrow(all_results), "个靶基因×细胞类型检验；",
    sum(all_results$nominal_p_lt_0_05), "项达到nominal P < 0.05。"
  ),
  paste0(
    "预设9基因中，",
    sum(all_results$prespecified_single_gene & all_results$nominal_p_lt_0_05),
    "项达到nominal P < 0.05。"
  ),
  paste0(
    sum(loo_summary$all_iterations_nominal_p_lt_0_05),
    "项主分析候选在逐一剔除16名ET供体的全部轮次中仍满足nominal P < 0.05。"
  ),
  "",
  "## 方法学依据",
  "",
  "- Castonguay et al. Nature Genetics. 2026. doi:10.1038/s41588-026-02544-8.",
  "- Law et al. Genome Biology. 2014. doi:10.1186/gb-2014-15-2-r29.",
  "- Ritchie et al. Nucleic Acids Research. 2015. doi:10.1093/nar/gkv007.",
  "- Crowell et al. Nature Communications. 2020. doi:10.1038/s41467-020-19894-4.",
  "- Squair et al. Nature Communications. 2021. doi:10.1038/s41467-021-25960-2."
)
writeLines(
  summary_lines,
  file.path(out_dir, "RAW_COUNT_PSEUDOBULK_SUMMARY_zh.md"),
  useBytes = TRUE
)
capture.output(sessionInfo(), file = file.path(out_dir, "session_info.txt"))
