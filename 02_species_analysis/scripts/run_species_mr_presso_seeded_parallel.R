#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1L || length(args) > 2L) {
  stop("Usage: Rscript run_species_mr_presso_seeded_parallel.R <run_root> [workers]")
}

root <- normalizePath(args[[1]], winslash = "/", mustWork = TRUE)
workers_requested <- if (length(args) == 2L) as.integer(args[[2]]) else 7L
if (!is.finite(workers_requested) || workers_requested < 1L) {
  stop("workers must be a positive integer")
}

results_dir <- file.path(root, "results")
logs_dir <- file.path(root, "logs")
input_path <- file.path(results_dir, "species_harmonised_instruments.csv")
output_path <- file.path(results_dir, "species_mr_presso_seeded.csv")
session_path <- file.path(results_dir, "r_session_info_seeded.txt")
determinism_path <- file.path(results_dir, "species_mr_presso_determinism_check.csv")

seed_base <- 20260902L
nb_distribution <- 10000L
group_columns <- c("gwas_accession", "taxon_name", "tier", "sensitivity_model")

derive_seed <- function(accession, tier, sensitivity_model) {
  key <- paste(accession, tier, sensitivity_model, sep = "|")
  code_points <- utf8ToInt(key)
  weighted_sum <- sum((seq_along(code_points) + 17) * code_points)
  as.integer((seed_base + weighted_sum) %% 2147483646 + 1)
}

run_one_task <- function(task) {
  dat <- task$data
  info <- task$info
  seed <- derive_seed(info$gwas_accession, info$tier, info$sensitivity_model)
  set.seed(
    seed,
    kind = "Mersenne-Twister",
    normal.kind = "Inversion",
    sample.kind = "Rejection"
  )
  rng_kind <- paste(RNGkind(), collapse = ";")

  result <- tryCatch(
    suppressMessages(MRPRESSO::mr_presso(
      BetaOutcome = "beta.outcome",
      BetaExposure = "beta.exposure",
      SdOutcome = "se.outcome",
      SdExposure = "se.exposure",
      OUTLIERtest = TRUE,
      DISTORTIONtest = TRUE,
      data = dat,
      NbDistribution = nb_distribution,
      SignifThreshold = 0.05
    )),
    error = identity
  )

  base <- data.frame(
    gwas_accession = info$gwas_accession,
    taxon_name = info$taxon_name,
    tier = info$tier,
    sensitivity_model = info$sensitivity_model,
    nsnp = nrow(dat),
    rng_seed = seed,
    rng_kind = rng_kind,
    nb_distribution = nb_distribution,
    stringsAsFactors = FALSE
  )

  if (inherits(result, "error")) {
    base$global_rss_observed <- NA_real_
    base$global_pvalue <- NA_character_
    base$outlier_count <- NA_integer_
    base$distortion_pvalue <- NA_character_
    base$analysis_status <- "not_estimable"
    base$error_message <- conditionMessage(result)
    return(base)
  }

  global <- result[["MR-PRESSO results"]][["Global Test"]]
  outlier <- result[["MR-PRESSO results"]][["Outlier Test"]]
  distortion <- result[["MR-PRESSO results"]][["Distortion Test"]]
  base$global_rss_observed <- global$RSSobs
  base$global_pvalue <- as.character(global$Pvalue)
  base$outlier_count <- if (is.null(outlier)) 0L else nrow(outlier)
  base$distortion_pvalue <- if (is.null(distortion)) {
    NA_character_
  } else {
    as.character(distortion$Pvalue)
  }
  base$analysis_status <- "estimated"
  base$error_message <- NA_character_
  base
}

if (!requireNamespace("MRPRESSO", quietly = TRUE)) {
  stop("MRPRESSO is not installed in the active R library")
}

dat <- read.csv(input_path, check.names = FALSE, stringsAsFactors = FALSE)
dat <- dat[dat$mr_keep %in% TRUE, , drop = FALSE]
keys <- unique(dat[group_columns])
keys <- keys[do.call(order, keys), , drop = FALSE]

tasks <- lapply(seq_len(nrow(keys)), function(index) {
  key <- keys[index, , drop = FALSE]
  keep <- rep(TRUE, nrow(dat))
  for (column in group_columns) {
    keep <- keep & dat[[column]] == key[[column]][1]
  }
  group_data <- dat[keep, , drop = FALSE]
  list(info = key, data = group_data)
})
tasks <- Filter(function(task) nrow(task$data) >= 4L, tasks)
if (length(tasks) == 0L) {
  stop("No harmonized species analysis had at least four instruments")
}

workers <- min(workers_requested, length(tasks))
cat("MR-PRESSO tasks:", length(tasks), " workers:", workers, "\n")
cluster <- parallel::makePSOCKcluster(
  workers,
  outfile = file.path(logs_dir, "mr_presso_worker_output.log")
)
on.exit({
  if (!is.null(cluster)) {
    parallel::stopCluster(cluster)
  }
}, add = TRUE)
parallel::clusterExport(
  cluster,
  c("derive_seed", "run_one_task", "seed_base", "nb_distribution"),
  envir = environment()
)
rows <- parallel::parLapplyLB(cluster, tasks, run_one_task)
parallel::stopCluster(cluster)
cluster <- NULL

table <- do.call(rbind, rows)
row.names(table) <- NULL
write.csv(table, output_path, row.names = FALSE)

check_row <- run_one_task(tasks[[1]])
reference_row <- table[1, names(check_row), drop = FALSE]
row.names(reference_row) <- NULL
row.names(check_row) <- NULL
check_pass <- isTRUE(all.equal(reference_row, check_row, check.attributes = FALSE))
write.csv(
  data.frame(
    gwas_accession = reference_row$gwas_accession,
    tier = reference_row$tier,
    sensitivity_model = reference_row$sensitivity_model,
    rng_seed = reference_row$rng_seed,
    exact_repeat_pass = check_pass,
    stringsAsFactors = FALSE
  ),
  determinism_path,
  row.names = FALSE
)
if (!check_pass) {
  stop("Determinism check failed for the first MR-PRESSO task")
}

sink(session_path)
sessionInfo()
sink()

cat("Seeded species MR-PRESSO completed:", output_path, "\n")
cat("Determinism check:", check_pass, "\n")
