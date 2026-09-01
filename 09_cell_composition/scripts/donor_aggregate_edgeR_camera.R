args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4) {
  stop(paste(
    "Usage: Rscript donor_aggregate_edgeR_camera.R",
    "<donor_count_dir> <target_matrix.csv> <sensitivity_out> <camera_out>"
  ))
}

suppressPackageStartupMessages(library(edgeR))
suppressPackageStartupMessages(library(limma))
options(warn = 1)

count_dir <- normalizePath(args[[1]], winslash = "/", mustWork = TRUE)
target_path <- normalizePath(args[[2]], winslash = "/", mustWork = TRUE)
sensitivity_out <- args[[3]]
camera_out <- args[[4]]
dir.create(sensitivity_out, recursive = TRUE, showWarnings = FALSE)
dir.create(camera_out, recursive = TRUE, showWarnings = FALSE)
target_model_dir <- file.path(sensitivity_out, "model_level_target_results")
ca_model_dir <- file.path(sensitivity_out, "model_level_carbonic_anhydrase_results")
camera_model_dir <- file.path(camera_out, "model_level_results")
dir.create(target_model_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(ca_model_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(camera_model_dir, recursive = TRUE, showWarnings = FALSE)

genes <- read.csv(file.path(count_dir, "genes.csv"), check.names = FALSE)
manifest <- read.csv(file.path(count_dir, "aggregation_manifest.csv"), check.names = FALSE)
targets <- read.csv(
  target_path,
  check.names = FALSE,
  fileEncoding = "UTF-8-BOM"
)
prediction_columns <- grep("_predicted$", names(targets), value = TRUE)
if (length(prediction_columns) != 8) stop("Expected eight chemical-form prediction columns")
targets[prediction_columns] <- lapply(
  targets[prediction_columns],
  function(values) tolower(as.character(values)) == "true"
)
conjugate_columns <- prediction_columns[
  grepl("SULFATE|GLUCURONIDE", prediction_columns)
]
if (length(conjugate_columns) != 4) {
  stop("Expected four sulfate or glucuronide prediction columns")
}

targets$predicted_for_any_chemical_form <-
  rowSums(targets[prediction_columns]) > 0
targets$predicted_for_all_chemical_forms <-
  rowSums(targets[prediction_columns]) == length(prediction_columns)
targets$predicted_for_any_sulfate_or_glucuronide <-
  rowSums(targets[conjugate_columns]) > 0
targets$predicted_for_all_sulfate_and_glucuronide_forms <-
  rowSums(targets[conjugate_columns]) == length(conjugate_columns)

target_flags <- targets[, c(
  "gene_symbol",
  "predicted_for_any_chemical_form",
  "predicted_for_all_chemical_forms",
  "predicted_for_any_sulfate_or_glucuronide",
  "predicted_for_all_sulfate_and_glucuronide_forms"
)]
if (nrow(target_flags) != 198 || anyDuplicated(target_flags$gene_symbol)) {
  stop("The prespecified target universe must contain 198 unique genes")
}

carbonic_anhydrase_family <- c(
  "CA1", "CA2", "CA3", "CA4", "CA5A", "CA5B", "CA6", "CA7",
  "CA8", "CA9", "CA10", "CA11", "CA12", "CA13", "CA14"
)
gene_sets <- list(
  "Targets predicted for any chemical form" = target_flags$gene_symbol[
    target_flags$predicted_for_any_chemical_form
  ],
  "Targets predicted for all chemical forms" = target_flags$gene_symbol[
    target_flags$predicted_for_all_chemical_forms
  ],
  "Targets predicted for any sulfate or glucuronide conjugate" =
    target_flags$gene_symbol[
      target_flags$predicted_for_any_sulfate_or_glucuronide
    ],
  "Targets predicted for all sulfate and glucuronide conjugates" =
    target_flags$gene_symbol[
      target_flags$predicted_for_all_sulfate_and_glucuronide_forms
    ],
  "Carbonic anhydrase family" = carbonic_anhydrase_family
)
expected_set_sizes <- c(198L, 25L, 150L, 45L, 15L)
if (!identical(as.integer(lengths(gene_sets)), expected_set_sizes)) {
  stop(paste(
    "Unexpected prespecified gene-set sizes:",
    paste(lengths(gene_sets), collapse = ",")
  ))
}

read_count_matrix <- function(slug, n_genes, n_donors) {
  path <- file.path(
    count_dir,
    "count_matrices",
    paste0(slug, ".gene_by_donor.int32.bin")
  )
  connection <- file(path, open = "rb")
  on.exit(close(connection))
  values <- readBin(
    connection,
    what = integer(),
    n = n_genes * n_donors,
    size = 4,
    signed = TRUE,
    endian = "little"
  )
  if (length(values) != n_genes * n_donors) {
    stop(paste("Unexpected binary length for", slug))
  }
  matrix(values, nrow = n_genes, ncol = n_donors)
}

prepare_design <- function(metadata) {
  metadata$condition <- factor(
    ifelse(metadata$disease == "essential tremor", "case", "control"),
    levels = c("control", "case")
  )
  metadata$age_c <- metadata$age - mean(metadata$age)
  metadata$sex_factor <- factor(metadata$sex, levels = c("female", "male"))
  metadata$batch_profile_factor <- factor(metadata$batch_profile)
  design <- model.matrix(
    ~ condition + age_c + sex_factor + batch_profile_factor,
    data = metadata
  )
  if (qr(design)$rank != ncol(design)) {
    stop("Donor-aggregate design matrix is not full rank")
  }
  list(metadata = metadata, design = design)
}

thresholds <- c(10L, 20L, 30L)
target_results <- list()
carbonic_anhydrase_results <- list()
camera_results <- list()
model_summaries <- list()
filter_audits <- list()

for (manifest_index in seq_len(nrow(manifest))) {
  cell_type <- manifest$cell_type[[manifest_index]]
  slug <- manifest$slug[[manifest_index]]
  metadata_all <- read.csv(
    file.path(count_dir, "sample_metadata", paste0(slug, ".csv")),
    check.names = FALSE,
    colClasses = c(donor_id = "character")
  )
  counts_all <- read_count_matrix(slug, nrow(genes), nrow(metadata_all))
  rownames(counts_all) <- genes$gene_symbol
  colnames(counts_all) <- metadata_all$donor_id
  if (any(counts_all < 0) || !all(colSums(counts_all) == metadata_all$library_size)) {
    stop(paste("Raw-count integrity check failed for", cell_type))
  }

  for (minimum_nuclei in thresholds) {
    keep_donors <- metadata_all$cell_count >= minimum_nuclei
    metadata <- metadata_all[keep_donors, , drop = FALSE]
    counts <- counts_all[, keep_donors, drop = FALSE]
    prepared <- prepare_design(metadata)
    metadata <- prepared$metadata
    design <- prepared$design

    dge <- DGEList(counts = counts)
    keep_genes <- filterByExpr(dge, design = design)
    if (sum(keep_genes) == 0) {
      stop(paste("No genes passed filterByExpr for", cell_type, minimum_nuclei))
    }
    dge <- dge[keep_genes, , keep.lib.sizes = FALSE]
    dge <- calcNormFactors(dge)
    dge <- estimateDisp(dge, design, robust = TRUE)
    ql_fit <- glmQLFit(dge, design, robust = TRUE)
    ql_test <- glmQLFTest(ql_fit, coef = "conditioncase")
    full_table <- topTags(
      ql_test,
      n = Inf,
      adjust.method = "BH",
      sort.by = "none"
    )$table
    full_table$gene_symbol <- rownames(full_table)
    rownames(full_table) <- NULL
    names(full_table)[names(full_table) == "FDR"] <-
      "BH_FDR_full_transcriptome_within_cell_type_supplementary"

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
    target_table$cell_type <- cell_type
    target_table$minimum_nuclei <- minimum_nuclei
    target_table$n_donors <- nrow(metadata)
    target_table$n_et_donors <- sum(metadata$condition == "case")
    target_table$n_control_donors <- sum(metadata$condition == "control")
    target_table$nominal_p_lt_0_05 <- target_table$PValue < 0.05
    target_results[[paste(cell_type, minimum_nuclei, sep = "|")]] <- target_table
    write.csv(
      target_table,
      file.path(target_model_dir, paste0(slug, "_", minimum_nuclei, "_nuclei.csv")),
      row.names = FALSE
    )

    ca_table <- full_table[
      full_table$gene_symbol %in% carbonic_anhydrase_family,
      ,
      drop = FALSE
    ]
    ca_table$cell_type <- cell_type
    ca_table$minimum_nuclei <- minimum_nuclei
    ca_table$n_donors <- nrow(metadata)
    ca_table$n_et_donors <- sum(metadata$condition == "case")
    ca_table$n_control_donors <- sum(metadata$condition == "control")
    ca_table$nominal_p_lt_0_05 <- ca_table$PValue < 0.05
    carbonic_anhydrase_results[[paste(cell_type, minimum_nuclei, sep = "|")]] <-
      ca_table
    write.csv(
      ca_table,
      file.path(ca_model_dir, paste0(slug, "_", minimum_nuclei, "_nuclei.csv")),
      row.names = FALSE
    )

    voom_object <- voom(dge, design, plot = FALSE)
    gene_set_indices <- ids2indices(
      gene_sets,
      rownames(voom_object$E),
      remove.empty = FALSE
    )
    for (correlation_method in c("fixed_0.01", "estimated")) {
      inter_gene_correlation <- if (
        correlation_method == "fixed_0.01"
      ) 0.01 else NA_real_
      camera_table <- camera(
        voom_object$E,
        index = gene_set_indices,
        design = design,
        contrast = "conditioncase",
        weights = voom_object$weights,
        inter.gene.cor = inter_gene_correlation,
        sort = FALSE
      )
      camera_table$gene_set <- rownames(camera_table)
      rownames(camera_table) <- NULL
      if (!"Correlation" %in% names(camera_table)) {
        camera_table$Correlation <- inter_gene_correlation
      }
      names(camera_table)[names(camera_table) == "FDR"] <-
        "BH_FDR_five_sets_within_cell_type_supplementary"
      camera_table$cell_type <- cell_type
      camera_table$minimum_nuclei <- minimum_nuclei
      camera_table$correlation_method <- correlation_method
      camera_table$nominal_p_lt_0_05 <- camera_table$PValue < 0.05
      camera_table$n_donors <- nrow(metadata)
      camera_table$n_et_donors <- sum(metadata$condition == "case")
      camera_table$n_control_donors <- sum(metadata$condition == "control")
      camera_results[[paste(
        cell_type,
        minimum_nuclei,
        correlation_method,
        sep = "|"
      )]] <- camera_table
      write.csv(
        camera_table,
        file.path(
          camera_model_dir,
          paste0(slug, "_", minimum_nuclei, "_nuclei_", correlation_method, ".csv")
        ),
        row.names = FALSE
      )
    }

    target_indices <- match(target_flags$gene_symbol, rownames(counts))
    filter_audits[[paste(cell_type, minimum_nuclei, sep = "|")]] <- data.frame(
      cell_type = cell_type,
      minimum_nuclei = minimum_nuclei,
      gene_symbol = target_flags$gene_symbol,
      total_raw_count = rowSums(counts[target_indices, , drop = FALSE]),
      donors_with_positive_count = rowSums(counts[target_indices, , drop = FALSE] > 0),
      passed_filterByExpr = keep_genes[target_indices]
    )
    model_summaries[[paste(cell_type, minimum_nuclei, sep = "|")]] <- data.frame(
      cell_type = cell_type,
      minimum_nuclei = minimum_nuclei,
      donors = nrow(metadata),
      et_donors = sum(metadata$condition == "case"),
      control_donors = sum(metadata$condition == "control"),
      genes_before_filter = nrow(counts),
      genes_after_filterByExpr = sum(keep_genes),
      design_rows = nrow(design),
      design_columns = ncol(design),
      design_rank = qr(design)$rank,
      minimum_library_size = min(colSums(counts)),
      maximum_library_size = max(colSums(counts))
    )
    message(
      "MODEL_COMPLETE cell_type=", cell_type,
      " minimum_nuclei=", minimum_nuclei,
      " donors=", nrow(metadata),
      " retained_genes=", sum(keep_genes)
    )
    rm(dge, ql_fit, ql_test, full_table, voom_object)
    invisible(gc(FALSE))
  }
  rm(counts_all)
  invisible(gc(FALSE))
}

target_results <- do.call(rbind, target_results)
rownames(target_results) <- NULL
target_results$BH_FDR_all_target_gene_cell_type_tests_within_threshold_supplementary <-
  ave(
    target_results$PValue,
    target_results$minimum_nuclei,
    FUN = function(values) p.adjust(values, method = "BH")
  )
target_results <- target_results[order(
  target_results$minimum_nuclei,
  target_results$PValue
), ]
write.csv(
  target_results,
  file.path(sensitivity_out, "all_target_gene_results_10_20_30_nuclei.csv"),
  row.names = FALSE
)
write.csv(
  target_results[target_results$nominal_p_lt_0_05, ],
  file.path(sensitivity_out, "nominal_target_candidates_10_20_30_nuclei.csv"),
  row.names = FALSE
)

carbonic_anhydrase_results <- do.call(rbind, carbonic_anhydrase_results)
rownames(carbonic_anhydrase_results) <- NULL
carbonic_anhydrase_results$BH_FDR_carbonic_anhydrase_tests_within_threshold_supplementary <-
  ave(
    carbonic_anhydrase_results$PValue,
    carbonic_anhydrase_results$minimum_nuclei,
    FUN = function(values) p.adjust(values, method = "BH")
  )
carbonic_anhydrase_results <- carbonic_anhydrase_results[order(
  carbonic_anhydrase_results$minimum_nuclei,
  carbonic_anhydrase_results$PValue
), ]
write.csv(
  carbonic_anhydrase_results,
  file.path(sensitivity_out, "carbonic_anhydrase_family_results_10_20_30_nuclei.csv"),
  row.names = FALSE
)

camera_results <- do.call(rbind, camera_results)
rownames(camera_results) <- NULL
camera_results$BH_FDR_70_tests_within_threshold_and_correlation_method_supplementary <-
  ave(
    camera_results$PValue,
    interaction(
      camera_results$minimum_nuclei,
      camera_results$correlation_method,
      drop = TRUE
    ),
    FUN = function(values) p.adjust(values, method = "BH")
  )
camera_results <- camera_results[order(
  camera_results$minimum_nuclei,
  camera_results$correlation_method,
  camera_results$PValue
), ]
write.csv(
  camera_results,
  file.path(camera_out, "camera_all_14_cell_types_10_20_30_nuclei.csv"),
  row.names = FALSE
)
write.csv(
  camera_results[camera_results$nominal_p_lt_0_05, ],
  file.path(camera_out, "camera_nominal_candidates.csv"),
  row.names = FALSE
)

camera_primary <- camera_results[
  camera_results$minimum_nuclei == 10 &
    camera_results$correlation_method == "fixed_0.01",
  ,
  drop = FALSE
]
write.csv(
  camera_primary,
  file.path(camera_out, "camera_primary_10_nuclei_all_70_tests.csv"),
  row.names = FALSE
)

model_summaries <- do.call(rbind, model_summaries)
rownames(model_summaries) <- NULL
write.csv(
  model_summaries,
  file.path(sensitivity_out, "model_summary_10_20_30_nuclei.csv"),
  row.names = FALSE
)
filter_audits <- do.call(rbind, filter_audits)
rownames(filter_audits) <- NULL
write.csv(
  filter_audits,
  file.path(sensitivity_out, "target_filter_audit_10_20_30_nuclei.csv"),
  row.names = FALSE
)

key_genes <- c("CA3", "CA4", "CA7", "CA12", "ADCY10", "TET3")
key_sensitivity <- target_results[
  target_results$gene_symbol %in% key_genes,
  c(
    "cell_type", "gene_symbol", "minimum_nuclei", "logFC", "PValue",
    "nominal_p_lt_0_05", "n_donors", "n_et_donors", "n_control_donors"
  ),
  drop = FALSE
]
write.csv(
  key_sensitivity,
  file.path(sensitivity_out, "key_gene_threshold_sensitivity.csv"),
  row.names = FALSE
)

summary_lines <- c(
  "# SCP3177供体聚合原始计数敏感性分析",
  "",
  "## 方法",
  "",
  "原始整数计数按供体×作者定义细胞类型聚合，并分别要求每个聚合至少含10、20或30个细胞核。",
  "每类细胞使用edgeR filterByExpr、TMM标准化和稳健quasi-likelihood模型；模型包含ET状态、年龄、性别和供体测序批次组合。",
  "候选纳入、保留和展示由nominal P < 0.05决定；BH校正值仅作为补充证据强度报告。",
  "CAMERA在相同供体聚合设计中完成5个预设集合×14类细胞的竞争性检验，并同时记录固定0.01及数据估计的基因间相关性设置。",
  "",
  "## 完成度",
  "",
  paste0("完成", nrow(model_summaries), "个细胞类型×最低细胞核阈值模型。"),
  paste0(
    "固定相关性设置下，10细胞核主CAMERA模型的70项检验中，",
    sum(camera_primary$nominal_p_lt_0_05), "项达到nominal P < 0.05。"
  ),
  "",
  "## 方法依据",
  "",
  "- Lun ATL, et al. Methods in Molecular Biology. 2016;1418:391-416. doi:10.1007/978-1-4939-3578-9_19.",
  "- Squair JW, et al. Nature Communications. 2021;12:5692. doi:10.1038/s41467-021-25960-2.",
  "- Wu D, Smyth GK. Nucleic Acids Research. 2012;40:e133. doi:10.1093/nar/gks461."
)
writeLines(
  summary_lines,
  file.path(sensitivity_out, "DONOR_AGGREGATE_RESULTS_zh.md"),
  useBytes = TRUE
)
capture.output(sessionInfo(), file = file.path(sensitivity_out, "session_info.txt"))
