suppressPackageStartupMessages({
  library(DESeq2)
  library(matrixStats)
})

command_arguments <- commandArgs(trailingOnly = FALSE)
script_argument <- grep("^--file=", command_arguments, value = TRUE)
if (length(script_argument) != 1L) {
  stop("Run this script with Rscript so that the data directory can be resolved.")
}
script_path <- normalizePath(sub("^--file=", "", script_argument), winslash = "/")
root <- dirname(dirname(script_path))
source_dir <- file.path(root, "09_source_data")
count_path <- file.path(source_dir, "gse197345_raw_count_matrix.csv.gz")
metadata_path <- file.path(source_dir, "gse197345_sample_metadata.csv")

counts_df <- read.csv(gzfile(count_path), check.names = FALSE)
rownames(counts_df) <- counts_df$gene_symbol
counts_df$gene_symbol <- NULL
counts <- as.matrix(counts_df)
storage.mode(counts) <- "integer"

metadata <- read.csv(metadata_path, stringsAsFactors = FALSE)
metadata$condition <- factor(metadata$condition, levels = c("Control", "ET"))
rownames(metadata) <- metadata$sample_id
if (!identical(colnames(counts), metadata$sample_id)) {
  stop("Count-matrix columns do not match sample metadata order.")
}

dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData = metadata,
  design = ~ condition
)
dds <- estimateSizeFactors(dds)
rld <- rlog(dds, blind = TRUE)
rlog_matrix <- assay(rld)

row_variance <- rowVars(rlog_matrix)
eligible <- which(is.finite(row_variance) & row_variance > 0)
top_n <- min(500L, length(eligible))
top_index <- eligible[order(row_variance[eligible], decreasing = TRUE)[seq_len(top_n)]]
pca <- prcomp(t(rlog_matrix[top_index, , drop = FALSE]), center = TRUE, scale. = FALSE)

# Principal-component signs are arbitrary. Each axis is oriented so that the
# loading with the largest absolute magnitude is positive.
for (axis in seq_len(min(2L, ncol(pca$rotation)))) {
  anchor <- which.max(abs(pca$rotation[, axis]))
  if (pca$rotation[anchor, axis] < 0) {
    pca$x[, axis] <- -pca$x[, axis]
    pca$rotation[, axis] <- -pca$rotation[, axis]
  }
}

variance_explained <- 100 * pca$sdev^2 / sum(pca$sdev^2)
pca_output <- data.frame(
  sample_id = rownames(pca$x),
  condition = metadata[rownames(pca$x), "condition"],
  PC1 = pca$x[, 1],
  PC2 = pca$x[, 2],
  stringsAsFactors = FALSE
)
write.csv(pca_output, file.path(source_dir, "gse197345_rlog_pca_coordinates.csv"), row.names = FALSE)
write.csv(
  data.frame(component = paste0("PC", seq_along(variance_explained)), variance_percent = variance_explained),
  file.path(source_dir, "gse197345_rlog_pca_variance.csv"),
  row.names = FALSE
)
write.csv(
  data.frame(gene_symbol = rownames(rlog_matrix)[top_index], variance = row_variance[top_index]),
  file.path(source_dir, "gse197345_rlog_pca_top500_genes.csv"),
  row.names = FALSE
)

heatmap_genes <- c("CA7", "ADCY10", "TET3", "KDM6B", "SIRT1", "IGF1R", "PTGES", "PTGS1")
missing_genes <- setdiff(heatmap_genes, rownames(rlog_matrix))
if (length(missing_genes) > 0L) {
  stop(paste("Heatmap genes missing:", paste(missing_genes, collapse = ", ")))
}
heatmap_matrix <- rlog_matrix[heatmap_genes, , drop = FALSE]
z_matrix <- t(scale(t(heatmap_matrix)))
if (any(!is.finite(z_matrix))) {
  stop("Non-finite z scores in the heatmap matrix.")
}

pca_by_sample <- pca_output
rownames(pca_by_sample) <- pca_by_sample$sample_id
sample_order <- unlist(lapply(levels(metadata$condition), function(group) {
  group_samples <- metadata$sample_id[metadata$condition == group]
  group_samples[order(pca_by_sample[group_samples, "PC1"])]
}), use.names = FALSE)

heatmap_output <- data.frame(
  gene_symbol = rownames(z_matrix),
  z_matrix[, sample_order, drop = FALSE],
  check.names = FALSE
)
write.csv(heatmap_output, file.path(source_dir, "gse197345_selected_gene_rlog_zscores.csv"), row.names = FALSE)
write.csv(
  data.frame(
    sample_id = sample_order,
    condition = metadata[sample_order, "condition"],
    display_order = seq_along(sample_order)
  ),
  file.path(source_dir, "gse197345_heatmap_sample_order.csv"),
  row.names = FALSE
)
write.csv(
  data.frame(sample_id = names(sizeFactors(dds)), size_factor = as.numeric(sizeFactors(dds))),
  file.path(source_dir, "gse197345_deseq2_size_factors.csv"),
  row.names = FALSE
)

cat(sprintf("Completed: %d genes, %d samples, %d PCA genes.\n", nrow(counts), ncol(counts), top_n))
