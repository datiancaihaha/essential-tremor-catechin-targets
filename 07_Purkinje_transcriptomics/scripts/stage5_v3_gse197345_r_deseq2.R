#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(DESeq2))

root <- "D:/CodexProjects/ET_MR_Stage0_20260821_v1"
tar_path <- file.path(
  root,
  "raw/stage5_cerebellum_validation_v2_20260822/GSE197345/GSE197345_RAW.tar"
)
work_dir <- file.path(root, "work/v3_parallel_evidence_20260822/GSE197345_counts")
out_dir <- file.path(root, "outputs/v3_parallel_evidence_20260822")
dir.create(work_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

untar(tar_path, exdir = work_dir)
files <- sort(list.files(work_dir, pattern = "\\.txt\\.gz$", full.names = TRUE))
stopifnot(length(files) == 40)

read_count <- function(path) {
  x <- read.delim(
    gzfile(path),
    header = FALSE,
    col.names = c("gene_symbol", "count"),
    colClasses = c("character", "integer"),
    check.names = FALSE
  )
  setNames(x$count, x$gene_symbol)
}

count_list <- lapply(files, read_count)
all_genes <- unique(unlist(lapply(count_list, names), use.names = FALSE))
count_matrix <- vapply(
  count_list,
  function(x) {
    y <- integer(length(all_genes))
    names(y) <- all_genes
    y[names(x)] <- x
    y
  },
  integer(length(all_genes))
)

sample_ids <- sub("_merged_markdups\\.txt\\.gz$", "", basename(files))
colnames(count_matrix) <- sample_ids
rownames(count_matrix) <- all_genes
count_matrix <- count_matrix[!startsWith(rownames(count_matrix), "__"), , drop = FALSE]
count_matrix <- count_matrix[rowSums(count_matrix) > 0, , drop = FALSE]

condition <- ifelse(grepl("_Control_", sample_ids), "Control", "ET")
col_data <- data.frame(
  condition = factor(condition, levels = c("Control", "ET")),
  row.names = sample_ids
)
stopifnot(sum(col_data$condition == "Control") == 16)
stopifnot(sum(col_data$condition == "ET") == 24)

dds <- DESeqDataSetFromMatrix(
  countData = count_matrix,
  colData = col_data,
  design = ~ condition
)
dds <- DESeq(dds)
res <- results(
  dds,
  contrast = c("condition", "ET", "Control"),
  alpha = 0.05
)

result <- as.data.frame(res)
result$gene_symbol <- rownames(result)
result$contrast <- "ET_vs_Control"
result$nominal_p_lt_0_05 <- !is.na(result$pvalue) & result$pvalue < 0.05
result$bh_fdr_lt_0_05_supplementary <- !is.na(result$padj) & result$padj < 0.05
result <- result[, c(
  "gene_symbol", "baseMean", "log2FoldChange", "lfcSE", "stat",
  "pvalue", "padj", "contrast", "nominal_p_lt_0_05",
  "bh_fdr_lt_0_05_supplementary"
)]
result <- result[order(result$pvalue, na.last = TRUE), ]

write.csv(
  result,
  file.path(out_dir, "08_gse197345_r_deseq2_all.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)
write.csv(
  data.frame(sample_id = rownames(col_data), condition = col_data$condition),
  file.path(out_dir, "09_gse197345_r_deseq2_sample_metadata.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

qa <- data.frame(
  metric = c(
    "r_version", "deseq2_version", "input_samples", "control_samples",
    "et_samples", "input_nonzero_biological_genes", "tested_genes",
    "nominal_p_lt_0_05", "bh_fdr_lt_0_05", "published_bh_fdr_lt_0_05"
  ),
  value = c(
    R.version.string,
    as.character(packageVersion("DESeq2")),
    ncol(count_matrix),
    sum(col_data$condition == "Control"),
    sum(col_data$condition == "ET"),
    nrow(count_matrix),
    sum(!is.na(result$pvalue)),
    sum(result$nominal_p_lt_0_05),
    sum(result$bh_fdr_lt_0_05_supplementary),
    36
  )
)
write.csv(
  qa,
  file.path(out_dir, "10_gse197345_r_deseq2_qa.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

message("R DESeq2 complete: ", nrow(result), " nonzero genes; ",
        sum(result$bh_fdr_lt_0_05_supplementary), " BH-FDR<0.05 genes")
