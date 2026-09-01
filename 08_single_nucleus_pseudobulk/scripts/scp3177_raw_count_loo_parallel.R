args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4) {
  stop(
    paste(
      "Usage: Rscript scp3177_raw_count_loo_parallel.R",
      "<pseudobulk_dir> <main_results.csv> <cell_type_slug|COMBINE> <out_dir>"
    )
  )
}

suppressPackageStartupMessages(library(edgeR))
suppressPackageStartupMessages(library(limma))
options(warn = 1)

pseudobulk_dir <- normalizePath(args[[1]], winslash = "/", mustWork = TRUE)
main_results_path <- normalizePath(args[[2]], winslash = "/", mustWork = TRUE)
requested_slug <- args[[3]]
out_dir <- args[[4]]
iteration_dir <- file.path(out_dir, "leave_one_et_donor_out_iterations")
cell_result_dir <- file.path(out_dir, "leave_one_et_donor_out_cell_type_results")
dir.create(iteration_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(cell_result_dir, recursive = TRUE, showWarnings = FALSE)

genes <- read.csv(file.path(pseudobulk_dir, "genes.csv"), check.names = FALSE)
manifest <- read.csv(
  file.path(pseudobulk_dir, "aggregation_manifest.csv"),
  check.names = FALSE
)
all_sample_metadata <- read.csv(
  file.path(pseudobulk_dir, "all_group_metadata.csv"),
  check.names = FALSE,
  colClasses = c(donor_id = "character", sample_id = "character")
)
main_results <- read.csv(main_results_path, check.names = FALSE)
seq_batch_levels <- sort(unique(all_sample_metadata$seq_batch))

if (requested_slug == "COMBINE") {
  expected_cell_types <- unique(main_results$cell_type[main_results$P.Value < 0.05])
  expected_slugs <- manifest$slug[match(expected_cell_types, manifest$cell_type)]
  result_paths <- file.path(
    cell_result_dir,
    paste0(expected_slugs, ".leave_one_et_donor_out.csv")
  )
  if (any(!file.exists(result_paths))) {
    stop(
      paste(
        "Missing cell-type leave-one-donor-out files:",
        paste(basename(result_paths[!file.exists(result_paths)]), collapse = ", ")
      )
    )
  }
  loo <- do.call(rbind, lapply(result_paths, function(path) {
    read.csv(path, check.names = FALSE, colClasses = c(omitted_et_donor = "character"))
  }))
  rownames(loo) <- NULL
  expected_rows <- sum(main_results$P.Value < 0.05) * 16L
  if (nrow(loo) != expected_rows) {
    stop(
      paste(
        "Unexpected combined leave-one-donor-out row count:",
        nrow(loo), "expected", expected_rows
      )
    )
  }
  key <- paste(loo$cell_type, loo$gene_symbol, loo$omitted_et_donor, sep = "|")
  if (anyDuplicated(key)) stop("Duplicate leave-one-donor-out result keys")
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
        maximum_nominal_p = if (all(is.na(values$P.Value))) {
          NA_real_
        } else {
          max(values$P.Value, na.rm = TRUE)
        },
        minimum_logFC = if (all(is.na(values$logFC))) {
          NA_real_
        } else {
          min(values$logFC, na.rm = TRUE)
        },
        maximum_logFC = if (all(is.na(values$logFC))) {
          NA_real_
        } else {
          max(values$logFC, na.rm = TRUE)
        }
      )
    }
  ))
  rownames(loo_summary) <- NULL
  write.csv(
    loo_summary,
    file.path(out_dir, "leave_one_et_donor_out_summary.csv"),
    row.names = FALSE
  )
  message(
    "COMBINE_COMPLETE rows=", nrow(loo),
    " candidates=", nrow(loo_summary),
    " all_iterations_nominal=",
    sum(loo_summary$all_iterations_nominal_p_lt_0_05)
  )
  quit(save = "no", status = 0)
}

manifest_row <- manifest[manifest$slug == requested_slug, , drop = FALSE]
if (nrow(manifest_row) != 1) stop("Requested cell-type slug is not unique in the manifest")
cell_type <- manifest_row$cell_type[[1]]
main_candidates <- main_results[
  main_results$cell_type == cell_type & main_results$P.Value < 0.05,
  c("gene_symbol", "logFC"),
  drop = FALSE
]
if (nrow(main_candidates) == 0) stop("Requested cell type has no nominal candidates")

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
    adjust.method = "BH"
  )
  list(table = table, keep = keep)
}

metadata <- read.csv(
  file.path(
    pseudobulk_dir,
    "sample_metadata",
    paste0(requested_slug, ".csv")
  ),
  check.names = FALSE,
  colClasses = c(donor_id = "character", sample_id = "character")
)
counts <- read_count_matrix(requested_slug, nrow(genes), nrow(metadata))
rownames(counts) <- genes$gene_symbol
colnames(counts) <- metadata$sample_id
if (!all(colSums(counts) == metadata$library_size)) {
  stop(paste("Library-size mismatch for", cell_type))
}
case_donors <- sort(unique(metadata$donor_id[metadata$disease == "essential tremor"]))
if (length(case_donors) != 16L) {
  stop(paste("Expected 16 ET donors for", cell_type, "but found", length(case_donors)))
}

iteration_paths <- character(length(case_donors))
for (donor_index in seq_along(case_donors)) {
  omitted_donor <- case_donors[[donor_index]]
  iteration_path <- file.path(
    iteration_dir,
    paste0(requested_slug, ".omit_", omitted_donor, ".csv")
  )
  iteration_paths[[donor_index]] <- iteration_path
  if (file.exists(iteration_path)) {
    prior <- read.csv(
      iteration_path,
      check.names = FALSE,
      colClasses = c(omitted_et_donor = "character")
    )
    if (
      nrow(prior) == nrow(main_candidates) &&
        identical(sort(prior$gene_symbol), sort(main_candidates$gene_symbol))
    ) {
      message(
        "LOO_RESUME cell_type=", cell_type,
        " donor=", donor_index, "/", length(case_donors)
      )
      next
    }
  }
  keep_samples <- metadata$donor_id != omitted_donor
  fitted <- fit_voom(
    counts[, keep_samples, drop = FALSE],
    metadata[keep_samples, , drop = FALSE]
  )
  table <- fitted$table
  table$gene_symbol <- rownames(table)
  matches <- match(main_candidates$gene_symbol, table$gene_symbol)
  iteration_result <- data.frame(
    cell_type = cell_type,
    gene_symbol = main_candidates$gene_symbol,
    omitted_et_donor = omitted_donor,
    passed_filterByExpr = !is.na(matches),
    logFC = ifelse(is.na(matches), NA_real_, table$logFC[matches]),
    P.Value = ifelse(is.na(matches), NA_real_, table$P.Value[matches]),
    main_logFC = main_candidates$logFC
  )
  iteration_result$same_direction_as_main <-
    sign(iteration_result$logFC) == sign(iteration_result$main_logFC)
  iteration_result$nominal_p_lt_0_05 <-
    !is.na(iteration_result$P.Value) & iteration_result$P.Value < 0.05
  write.csv(iteration_result, iteration_path, row.names = FALSE)
  message(
    "LOO_COMPLETE cell_type=", cell_type,
    " donor=", donor_index, "/", length(case_donors)
  )
  rm(fitted, table, iteration_result)
  invisible(gc(FALSE))
}

cell_results <- do.call(rbind, lapply(iteration_paths, function(path) {
  read.csv(path, check.names = FALSE, colClasses = c(omitted_et_donor = "character"))
}))
rownames(cell_results) <- NULL
write.csv(
  cell_results,
  file.path(
    cell_result_dir,
    paste0(requested_slug, ".leave_one_et_donor_out.csv")
  ),
  row.names = FALSE
)
message(
  "CELL_TYPE_LOO_COMPLETE cell_type=", cell_type,
  " rows=", nrow(cell_results)
)
