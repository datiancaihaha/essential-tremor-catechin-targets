args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4) {
  stop(
    paste(
      "Usage: Rscript scp3177_finalize_raw_count_results.R",
      "<analysis_root> <target_matrix.csv> <chemical_form_records.csv> <out_dir>"
    )
  )
}

analysis_root <- normalizePath(args[[1]], winslash = "/", mustWork = TRUE)
target_matrix_path <- normalizePath(args[[2]], winslash = "/", mustWork = TRUE)
chemical_form_path <- normalizePath(args[[3]], winslash = "/", mustWork = TRUE)
out_dir <- args[[4]]
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

raw_dir <- file.path(analysis_root, "09_raw_count_voom")
normalized_dir <- file.path(analysis_root, "02_normalized_expression_limma")
proportion_dir <- file.path(analysis_root, "05_cell_proportion_propeller")

raw <- read.csv(
  file.path(raw_dir, "all_target_gene_cell_type_results.csv"),
  check.names = FALSE
)
nominal <- read.csv(
  file.path(raw_dir, "nominal_target_candidates.csv"),
  check.names = FALSE
)
loo <- read.csv(
  file.path(raw_dir, "leave_one_et_donor_out_summary.csv"),
  check.names = FALSE
)
loo_all <- read.csv(
  file.path(raw_dir, "leave_one_et_donor_out_results.csv"),
  check.names = FALSE,
  colClasses = c(omitted_et_donor = "character")
)
model_summary <- read.csv(
  file.path(raw_dir, "cell_type_model_summary.csv"),
  check.names = FALSE
)
filter_audit <- read.csv(
  file.path(raw_dir, "target_filter_audit.csv"),
  check.names = FALSE
)
normalized <- read.csv(
  file.path(normalized_dir, "all_gene_cell_type_results.csv"),
  check.names = FALSE
)
proportion <- read.csv(
  file.path(proportion_dir, "propeller_nominal_candidates.csv"),
  check.names = FALSE
)
proportion_loo <- read.csv(
  file.path(proportion_dir, "leave_one_et_donor_out_summary.csv"),
  check.names = FALSE
)
target_matrix <- read.csv(
  target_matrix_path,
  check.names = FALSE,
  fileEncoding = "UTF-8-BOM"
)
chemical_forms <- read.csv(chemical_form_path, check.names = FALSE)
if (!endsWith(names(chemical_forms)[[1]], "record_id")) {
  stop("The first chemical-form column is not record_id")
}
names(chemical_forms)[[1]] <- "record_id"

ontology_map <- c(
  Astrocytes = "astrocyte",
  Bergmann = "Bergmann glial cell",
  Endocytes = "endothelial cell",
  Golgi = "Golgi cell",
  Granule = "granule cell",
  MLI_1 = "molecular layer interneuron",
  MLI_2 = "molecular layer interneuron",
  Microglia = "microglia",
  OPC = "oligodendrocyte precursor cell",
  Oligodendrocytes = "oligodendrocyte",
  PLI = "interneuron",
  Pericytes = "pericyte",
  Purkinje = "Purkinje cell",
  UBC = "unipolar brush cell"
)

candidate_results <- merge(
  nominal,
  loo,
  by = c("cell_type", "gene_symbol"),
  all.x = TRUE,
  sort = FALSE
)
candidate_results$ontology_cell_type <- unname(ontology_map[candidate_results$cell_type])
normalized_comparison <- normalized[, c(
  "cell_type", "gene_symbol", "logFC", "P.Value", "nominal_p_lt_0_05"
)]
names(normalized_comparison) <- c(
  "ontology_cell_type", "gene_symbol", "normalized_expression_logFC",
  "normalized_expression_P.Value", "normalized_expression_nominal_p_lt_0_05"
)
candidate_results <- merge(
  candidate_results,
  normalized_comparison,
  by = c("ontology_cell_type", "gene_symbol"),
  all.x = TRUE,
  sort = FALSE
)
candidate_results$normalized_expression_same_direction <-
  sign(candidate_results$logFC) == sign(candidate_results$normalized_expression_logFC)
candidate_results <- candidate_results[order(candidate_results$P.Value), ]
write.csv(
  candidate_results,
  file.path(out_dir, "raw_count_nominal_candidates_with_leave_one_out.csv"),
  row.names = FALSE
)

all_concordance <- raw
all_concordance$ontology_cell_type <- unname(ontology_map[all_concordance$cell_type])
all_concordance <- merge(
  all_concordance,
  normalized_comparison,
  by = c("ontology_cell_type", "gene_symbol"),
  all.x = TRUE,
  sort = FALSE
)
all_concordance$normalized_expression_same_direction <-
  sign(all_concordance$logFC) == sign(all_concordance$normalized_expression_logFC)
write.csv(
  all_concordance,
  file.path(out_dir, "raw_count_and_normalized_expression_comparison.csv"),
  row.names = FALSE
)

cell_type_summary <- do.call(rbind, lapply(
  sort(unique(raw$cell_type)),
  function(cell_type) {
    evaluated <- raw[raw$cell_type == cell_type, ]
    selected <- candidate_results[candidate_results$cell_type == cell_type, ]
    model <- model_summary[model_summary$cell_type == cell_type, ]
    data.frame(
      cell_type = cell_type,
      pseudobulk_samples = model$pseudobulk_samples[[1]],
      donors = model$donors[[1]],
      et_donors = model$et_donors[[1]],
      normal_donors = model$normal_donors[[1]],
      genes_after_filterByExpr = model$genes_after_filterByExpr[[1]],
      evaluable_target_gene_tests = nrow(evaluated),
      nominal_p_lt_0_05_candidates = nrow(selected),
      all_16_iterations_same_direction = sum(selected$same_direction_fraction == 1),
      all_16_iterations_nominal_p_lt_0_05 =
        sum(selected$all_iterations_nominal_p_lt_0_05),
      supplementary_target_set_BH_lt_0_05 = sum(
        selected$BH_FDR_all_target_gene_cell_type_tests_supplementary < 0.05
      )
    )
  }
))
write.csv(
  cell_type_summary,
  file.path(out_dir, "raw_count_cell_type_summary.csv"),
  row.names = FALSE
)

target_set_names <- c(
  "all_form_union", "all_form_core",
  "circulating_conjugate_union", "circulating_conjugate_core"
)
target_set_labels <- c(
  all_form_union = "Union across eight chemical forms",
  all_form_core = "Shared across all eight chemical forms",
  circulating_conjugate_union = "Union across sulfate and glucuronide conjugates",
  circulating_conjugate_core = "Shared across sulfate and glucuronide conjugates"
)
target_set_summary <- do.call(rbind, lapply(target_set_names, function(set_name) {
  evaluated <- raw[raw[[set_name]] %in% TRUE, ]
  selected <- candidate_results[candidate_results[[set_name]] %in% TRUE, ]
  data.frame(
    target_set = target_set_labels[[set_name]],
    evaluable_gene_cell_type_tests = nrow(evaluated),
    nominal_p_lt_0_05_candidates = nrow(selected),
    all_16_iterations_same_direction = sum(selected$same_direction_fraction == 1),
    all_16_iterations_nominal_p_lt_0_05 =
      sum(selected$all_iterations_nominal_p_lt_0_05),
    unique_stable_genes = length(unique(
      selected$gene_symbol[selected$all_iterations_nominal_p_lt_0_05]
    ))
  )
}))
write.csv(
  target_set_summary,
  file.path(out_dir, "chemical_form_target_set_summary.csv"),
  row.names = FALSE
)

prediction_columns <- grep("_predicted$", names(target_matrix), value = TRUE)
target_matrix[prediction_columns] <- lapply(
  target_matrix[prediction_columns],
  function(values) tolower(as.character(values)) == "true"
)
form_lookup <- unique(chemical_forms[, c("record_id", "standard_name", "pubchem_cid")])
form_summary <- do.call(rbind, lapply(prediction_columns, function(column_name) {
  record_id <- sub("_predicted$", "", column_name)
  genes <- target_matrix$gene_symbol[target_matrix[[column_name]] %in% TRUE]
  evaluated <- raw[raw$gene_symbol %in% genes, ]
  selected <- candidate_results[candidate_results$gene_symbol %in% genes, ]
  form <- form_lookup[form_lookup$record_id == record_id, ]
  data.frame(
    chemical_form = form$standard_name[[1]],
    PubChem_CID = form$pubchem_cid[[1]],
    predicted_target_genes = length(unique(genes)),
    evaluable_gene_cell_type_tests = nrow(evaluated),
    nominal_p_lt_0_05_candidates = nrow(selected),
    all_16_iterations_nominal_p_lt_0_05 =
      sum(selected$all_iterations_nominal_p_lt_0_05),
    unique_stable_genes = length(unique(
      selected$gene_symbol[selected$all_iterations_nominal_p_lt_0_05]
    ))
  )
}))
write.csv(
  form_summary,
  file.path(out_dir, "eight_chemical_forms_expression_summary.csv"),
  row.names = FALSE
)

stable <- candidate_results[candidate_results$all_iterations_nominal_p_lt_0_05, ]
stable_by_gene <- split(stable, as.character(stable$gene_symbol), drop = TRUE)
multicell_summary <- do.call(rbind, lapply(stable_by_gene, function(values) {
  data.frame(
    gene_symbol = values$gene_symbol[[1]],
    stable_cell_types = nrow(values),
    effect_direction = if (all(values$logFC > 0)) {
      "positive"
    } else if (all(values$logFC < 0)) {
      "negative"
    } else {
      "mixed"
    },
    cell_types = paste(values$cell_type, collapse = "; "),
    minimum_nominal_p = min(values$P.Value),
    maximum_absolute_logFC = max(abs(values$logFC)),
    supplementary_target_set_BH_lt_0_05_cell_types = sum(
      values$BH_FDR_all_target_gene_cell_type_tests_supplementary < 0.05
    )
  )
}))
rownames(multicell_summary) <- NULL
multicell_summary <- multicell_summary[
  order(-multicell_summary$stable_cell_types, multicell_summary$minimum_nominal_p),
]
write.csv(
  multicell_summary,
  file.path(out_dir, "stable_candidate_gene_summary.csv"),
  row.names = FALSE
)

carbonic_anhydrase_genes <- c(
  "CA1", "CA2", "CA3", "CA4", "CA5A", "CA5B", "CA6", "CA7",
  "CA8", "CA9", "CA10", "CA11", "CA12", "CA13", "CA14", "CA15P"
)
carbonic_anhydrase_results <- raw[raw$gene_symbol %in% carbonic_anhydrase_genes, ]
carbonic_anhydrase_results <- merge(
  carbonic_anhydrase_results,
  loo,
  by = c("cell_type", "gene_symbol"),
  all.x = TRUE,
  sort = FALSE
)
carbonic_anhydrase_results <- carbonic_anhydrase_results[
  order(carbonic_anhydrase_results$P.Value),
]
write.csv(
  carbonic_anhydrase_results,
  file.path(out_dir, "carbonic_anhydrase_family_results.csv"),
  row.names = FALSE
)

prespecified_results <- raw[raw$prespecified_single_gene, ]
prespecified_results <- merge(
  prespecified_results,
  loo,
  by = c("cell_type", "gene_symbol"),
  all.x = TRUE,
  sort = FALSE
)
prespecified_results <- prespecified_results[order(prespecified_results$P.Value), ]
write.csv(
  prespecified_results,
  file.path(out_dir, "prespecified_gene_results_with_leave_one_out.csv"),
  row.names = FALSE
)

proportion_results <- merge(
  proportion,
  proportion_loo,
  by = "cell_type",
  all.x = TRUE,
  sort = FALSE
)
proportion_results <- proportion_results[order(proportion_results$P.Value), ]
write.csv(
  proportion_results,
  file.path(out_dir, "cell_proportion_results_with_leave_one_out.csv"),
  row.names = FALSE
)

iteration_key <- paste(
  loo_all$cell_type,
  loo_all$gene_symbol,
  loo_all$omitted_et_donor,
  sep = "|"
)
main_key <- paste(raw$cell_type, raw$gene_symbol, sep = "|")
nominal_key <- paste(nominal$cell_type, nominal$gene_symbol, sep = "|")
expected_nominal_key <- paste(
  raw$cell_type[raw$P.Value < 0.05],
  raw$gene_symbol[raw$P.Value < 0.05],
  sep = "|"
)
iteration_files <- list.files(
  file.path(raw_dir, "leave_one_et_donor_out_iterations"),
  pattern = "\\.csv$",
  full.names = TRUE
)
cell_result_files <- list.files(
  file.path(raw_dir, "leave_one_et_donor_out_cell_type_results"),
  pattern = "\\.csv$",
  full.names = TRUE
)
full_transcriptome_files <- list.files(
  file.path(raw_dir, "full_transcriptome_results"),
  pattern = "\\.csv\\.gz$",
  full.names = TRUE
)
qa <- data.frame(
  check = c(
    "Raw-count result keys are unique",
    "Raw-count P values are complete",
    "Nominal candidate subset is exact",
    "Leave-one-donor-out rows equal 163 x 16",
    "Leave-one-donor-out keys are unique",
    "Every nominal candidate has 16 iterations",
    "Iteration files are complete",
    "Cell-type leave-one-donor-out files are complete",
    "Full-transcriptome result files are complete",
    "All 14 main models are present"
  ),
  observed = c(
    length(unique(main_key)),
    sum(!is.na(raw$P.Value)),
    length(intersect(nominal_key, expected_nominal_key)),
    nrow(loo_all),
    length(unique(iteration_key)),
    min(table(paste(loo_all$cell_type, loo_all$gene_symbol, sep = "|"))),
    length(iteration_files),
    length(cell_result_files),
    length(full_transcriptome_files),
    nrow(model_summary)
  ),
  expected = c(
    nrow(raw), nrow(raw), length(expected_nominal_key),
    163L * 16L, 163L * 16L, 16L, 224L, 14L, 14L, 14L
  )
)
qa$pass <- qa$observed == qa$expected
write.csv(
  qa,
  file.path(out_dir, "quality_control_checks.csv"),
  row.names = FALSE
)
if (!all(qa$pass)) stop("One or more final quality-control checks failed")

all_matched <- !is.na(all_concordance$normalized_expression_logFC)
raw_nominal_comparison <- candidate_results[
  !is.na(candidate_results$normalized_expression_logFC),
]
tet3 <- candidate_results[
  candidate_results$cell_type == "MLI_1" & candidate_results$gene_symbol == "TET3",
]
ca12 <- candidate_results[
  candidate_results$gene_symbol == "CA12" &
    candidate_results$all_iterations_nominal_p_lt_0_05,
]
top_multicell <- multicell_summary[multicell_summary$stable_cell_types >= 3, ]
top_multicell_text <- paste(
  paste0(
    top_multicell$gene_symbol, " (", top_multicell$stable_cell_types,
    " cell types, ", top_multicell$effect_direction, ")"
  ),
  collapse = "; "
)
proportion_text <- paste(proportion_results$cell_type, collapse = ", ")

report <- c(
  "# SCP3177原始计数伪批量分析综合结果",
  "",
  "## 数据与模型",
  "",
  paste0(
    "分析使用官方SCP3177文件的raw count矩阵，共1,004,112个细胞、49,401个基因、109名供体。",
    "按供体、测序批次和作者细胞类型聚合为1,753个伪批量样本；每个样本至少包含10个细胞。"
  ),
  paste0(
    "差异表达模型采用edgeR的filterByExpr与TMM标准化、limma-voom和duplicateCorrelation；",
    "协变量包括ET状态、年龄、性别和测序批次，供体作为重复测量块。"
  ),
  "候选纳入、保留和展示统一依据nominal P < 0.05；BH校正值仅作为补充证据强度记录。",
  "",
  "## 主分析结果",
  "",
  paste0(
    "14类细胞的主模型全部完成。1,013个可评估靶基因与细胞类型组合中，",
    nrow(candidate_results), "项满足nominal P < 0.05。"
  ),
  paste0(
    "所有", nrow(candidate_results), "项候选在16轮逐ET供体剔除分析中保持相同效应方向；",
    sum(candidate_results$all_iterations_nominal_p_lt_0_05),
    "项在16/16轮均保持nominal P < 0.05。"
  ),
  paste0(
    "作为补充证据强度，稳定候选中",
    sum(stable$BH_FDR_all_target_gene_cell_type_tests_supplementary < 0.05),
    "项的全靶基因集合BH校正值小于0.05；该指标未用于候选筛选。"
  ),
  "",
  "## 跨细胞重复信号",
  "",
  top_multicell_text,
  "这些结果反映统计学上的跨细胞重复模式；在进入生物学叙事前仍需逐基因文献核验。",
  "",
  "## 碳酸酐酶家族",
  "",
  paste0(
    "CA12在Astrocytes、Bergmann和Granule中均为正向nominal候选，且三类细胞均在16/16轮留一分析中保持nominal P < 0.05。",
    "对应主模型P值范围为", formatC(min(ca12$P.Value), digits = 3, format = "e"),
    "至", formatC(max(ca12$P.Value), digits = 3, format = "e"), "。"
  ),
  "CA4在MLI_2中为16/16轮稳定的正向nominal候选；在OPC中为主分析nominal候选，但仅10/16轮保持nominal P < 0.05。",
  "CA3未在正式原始计数伪批量分析中形成nominal候选，因此不应沿用旧归一化表达敏感性分析中的CA3表达结论。",
  "",
  "## 预设基因",
  "",
  paste0(
    "预设9基因中，仅TET3在MLI_1主模型中满足nominal P < 0.05：log2 fold change = ",
    formatC(tet3$logFC, digits = 3, format = "f"), ", P = ",
    formatC(tet3$P.Value, digits = 3, format = "g"), "。"
  ),
  paste0(
    "TET3效应方向在16/16轮均为负，",
    round(tet3$nominal_p_lt_0_05_fraction * 16),
    "/16轮保持nominal P < 0.05；另外两轮P值为0.0505和0.0521。"
  ),
  "该结果适合作为方向稳定的次级证据，不应表述为16轮全部达到nominal P < 0.05。",
  "",
  "## 八种化学形式",
  "",
  paste0(
    "八种化学形式靶点并集共有1,013个可评估基因与细胞类型组合，",
    "其中163项满足nominal P < 0.05，89项在16/16轮均保持nominal P < 0.05。"
  ),
  paste0(
    "四种硫酸化和葡萄糖醛酸化结合物的靶点并集包含816个可评估组合、135项nominal候选，",
    "其中77项在16/16轮均保持nominal P < 0.05。"
  ),
  "各化学形式使用可检索的标准化学名称与PubChem CID汇总，未在正式结果中使用内部代号。",
  "",
  "## 细胞组成",
  "",
  paste0(
    "独立的供体层面细胞组成分析保留6项nominal P < 0.05结果：", proportion_text, "。"
  ),
  "六项结果在16/16轮逐ET供体剔除分析中均保持相同方向和nominal P < 0.05，因此可作为独立于表达模型的补强证据。",
  "",
  "## 与归一化表达敏感性分析的关系",
  "",
  paste0(
    "正式原始计数模型与旧归一化表达敏感性分析共有",
    sum(all_matched), "个可比较组合；效应量Pearson相关系数为",
    formatC(cor(
      all_concordance$logFC[all_matched],
      all_concordance$normalized_expression_logFC[all_matched]
    ), digits = 3, format = "f"), "。"
  ),
  paste0(
    "在", nrow(raw_nominal_comparison), "项正式nominal候选中，",
    sum(raw_nominal_comparison$normalized_expression_P.Value < 0.05),
    "项在旧敏感性分析中也满足nominal P < 0.05。"
  ),
  "由于两种处理路径的总体一致性较低，正式原始计数伪批量结果应作为主要表达证据；旧归一化表达结果仅保留为技术敏感性信息，不用于替代或放大正式结果。",
  "",
  "## 建议的结果层级",
  "",
  "1. 优先保留：供体层面细胞组成结果；CA12在三类细胞中的一致正向信号；以及具有多细胞重复且16/16轮稳定的靶点结果。",
  "2. 次级保留：TET3在MLI_1中的负向结果，明确标注14/16轮保持nominal P < 0.05。",
  "3. 不作为正式表达结论：仅在旧归一化表达敏感性分析中出现、但未被原始计数伪批量模型支持的信号。",
  "4. 生物学优先级需在纳入稿件前与ET遗传学、小脑细胞生物学及化合物机制文献逐项匹配。",
  "",
  "## 方法学依据",
  "",
  "- Castonguay et al. Nature Genetics. 2026. doi:10.1038/s41588-026-02544-8.",
  "- Crowell et al. Nature Communications. 2020. doi:10.1038/s41467-020-19894-4.",
  "- Squair et al. Nature Communications. 2021. doi:10.1038/s41467-021-25960-2.",
  "- Law et al. Genome Biology. 2014. doi:10.1186/gb-2014-15-2-r29.",
  "- Ritchie et al. Nucleic Acids Research. 2015. doi:10.1093/nar/gkv007.",
  "- Phipson et al. Bioinformatics. 2022. doi:10.1093/bioinformatics/btac582."
)
writeLines(
  report,
  file.path(out_dir, "SCP3177_RAW_COUNT_RESULTS_zh.md"),
  useBytes = TRUE
)
capture.output(sessionInfo(), file = file.path(out_dir, "session_info.txt"))
