#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(TwoSampleMR)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) {
  stop("Usage: Rscript run_species_mr_seeded.R <run_root>")
}
root <- normalizePath(args[[1]], winslash = "/", mustWork = TRUE)
results_dir <- file.path(root, "results")
input_path <- file.path(results_dir, "species_mr_harmonization_input.csv")
crosswalk_path <- file.path(results_dir, "microbiome_taxonomic_crosswalk.csv")
presso_output_path <- file.path(results_dir, "species_mr_presso.csv")
run_presso <- TRUE
mr_presso_seed_base <- 20260902L
mr_presso_nb_distribution <- 10000L

et_prevalence <- 0.0032
et_cases <- 16480
et_controls <- 1936173
exposure_n <- 16017

tiers <- list(
  study_wide = "included_study_wide",
  genome_wide = "included_genome_wide",
  exploratory = "included_exploratory"
)

as_flag <- function(x) {
  tolower(as.character(x)) %in% c("true", "t", "1")
}

append_row <- function(store, row) {
  store[[length(store) + 1L]] <- row
  store
}

bind_rows_base <- function(rows) {
  if (length(rows) == 0L) {
    return(data.frame())
  }
  columns <- unique(unlist(lapply(rows, names), use.names = FALSE))
  normalized <- lapply(rows, function(row) {
    missing <- setdiff(columns, names(row))
    for (column in missing) {
      row[[column]] <- NA
    }
    row[columns]
  })
  do.call(rbind, normalized)
}

format_exposure <- function(x, tier, sensitivity_model) {
  prevalence <- unique(x$prevalence)
  prevalence <- prevalence[is.finite(prevalence)][1]
  ncase <- round(prevalence * exposure_n)
  data.frame(
    SNP = x$rs_id,
    beta.exposure = x$beta,
    se.exposure = x$standard_error,
    effect_allele.exposure = toupper(x$effect_allele),
    other_allele.exposure = toupper(x$other_allele),
    eaf.exposure = x$effect_allele_frequency,
    pval.exposure = x$p_value,
    samplesize.exposure = exposure_n,
    ncase.exposure = ncase,
    ncontrol.exposure = exposure_n - ncase,
    prevalence.exposure = prevalence,
    units.exposure = "log odds",
    exposure = unique(x$taxon_name)[1],
    id.exposure = paste(unique(x$gwas_accession)[1], tier, sensitivity_model, sep = "__"),
    stringsAsFactors = FALSE
  )
}

format_outcome <- function(x) {
  data.frame(
    SNP = x$rs_id,
    beta.outcome = x$outcome_beta,
    se.outcome = x$outcome_se,
    effect_allele.outcome = toupper(x$outcome_effect_allele),
    other_allele.outcome = toupper(x$outcome_other_allele),
    eaf.outcome = x$outcome_eaf_median,
    pval.outcome = x$outcome_p,
    samplesize.outcome = et_cases + et_controls,
    ncase.outcome = et_cases,
    ncontrol.outcome = et_controls,
    prevalence.outcome = et_prevalence,
    units.outcome = "log odds",
    outcome = "Essential tremor",
    id.outcome = "deCODE_essential_tremor",
    stringsAsFactors = FALSE
  )
}

safe_analysis <- function(expression) {
  tryCatch(
    list(status = "estimated", value = force(expression), message = NA_character_),
    error = function(error) {
      list(status = "not_estimable", value = NULL, message = conditionMessage(error))
    }
  )
}

derive_mr_presso_seed <- function(accession, tier, sensitivity_model) {
  key <- paste(accession, tier, sensitivity_model, sep = "|")
  code_points <- utf8ToInt(key)
  weighted_sum <- sum((seq_along(code_points) + 17) * code_points)
  as.integer((mr_presso_seed_base + weighted_sum) %% 2147483646 + 1)
}

input <- read.csv(input_path, check.names = FALSE, stringsAsFactors = FALSE)
crosswalk <- read.csv(crosswalk_path, check.names = FALSE, stringsAsFactors = FALSE)
input <- input[input$outcome_selection_status == "standard_allele_match", , drop = FALSE]

for (column in unlist(tiers, use.names = FALSE)) {
  input[[column]] <- as_flag(input[[column]])
}

analysis_status <- list()
harmonised_rows <- list()
mr_rows <- list()
heterogeneity_rows <- list()
pleiotropy_rows <- list()
steiger_rows <- list()
leave_one_out_rows <- list()
presso_rows <- list()

accessions <- unique(crosswalk$accession)

for (accession in accessions) {
  accession_data <- input[input$gwas_accession == accession, , drop = FALSE]
  species_name <- crosswalk$taxon_name[crosswalk$accession == accession][1]
  for (tier_name in names(tiers)) {
    tier_column <- tiers[[tier_name]]
    tier_data <- accession_data[accession_data[[tier_column]], , drop = FALSE]
    models <- list(full = tier_data)
    if (nrow(tier_data) > 0L && any(tier_data$explicit_host_locus == "ABO", na.rm = TRUE)) {
      models$exclude_ABO <- tier_data[tier_data$explicit_host_locus != "ABO", , drop = FALSE]
    }
    if (nrow(tier_data) > 0L && any(tier_data$explicit_host_locus == "FOXP1", na.rm = TRUE)) {
      models$exclude_FOXP1 <- tier_data[tier_data$explicit_host_locus != "FOXP1", , drop = FALSE]
    }
    if (nrow(tier_data) > 0L && any(
      tier_data$explicit_host_locus %in% c("ABO", "FUT2", "LCT", "FOXP1"),
      na.rm = TRUE
    )) {
      models$exclude_host_loci <- tier_data[
        !tier_data$explicit_host_locus %in% c("ABO", "FUT2", "LCT", "FOXP1"),
        ,
        drop = FALSE
      ]
    }

    for (model_name in names(models)) {
      model_data <- models[[model_name]]
      status_base <- data.frame(
        gwas_accession = accession,
        taxon_name = species_name,
        tier = tier_name,
        sensitivity_model = model_name,
        instruments_before_harmonisation = nrow(model_data),
        instruments_after_harmonisation = 0L,
        analysis_status = "no_instruments",
        message = NA_character_,
        stringsAsFactors = FALSE
      )
      if (nrow(model_data) == 0L) {
        analysis_status <- append_row(analysis_status, status_base)
        next
      }

      exposure <- format_exposure(model_data, tier_name, model_name)
      outcome <- format_outcome(model_data)
      harmonised_result <- safe_analysis(
        suppressMessages(harmonise_data(exposure, outcome, action = 2))
      )
      if (harmonised_result$status != "estimated") {
        status_base$analysis_status <- harmonised_result$status
        status_base$message <- harmonised_result$message
        analysis_status <- append_row(analysis_status, status_base)
        next
      }
      dat <- harmonised_result$value
      dat <- dat[dat$mr_keep %in% TRUE, , drop = FALSE]
      status_base$instruments_after_harmonisation <- nrow(dat)
      if (nrow(dat) == 0L) {
        status_base$analysis_status <- "no_harmonised_instruments"
        analysis_status <- append_row(analysis_status, status_base)
        next
      }
      status_base$analysis_status <- "estimated"
      analysis_status <- append_row(analysis_status, status_base)

      dat$gwas_accession <- accession
      dat$taxon_name <- species_name
      dat$tier <- tier_name
      dat$sensitivity_model <- model_name
      dat$explicit_host_locus <- model_data$explicit_host_locus[match(dat$SNP, model_data$rs_id)]
      harmonised_rows <- append_row(harmonised_rows, dat)

      methods <- if (nrow(dat) == 1L) {
        "mr_wald_ratio"
      } else if (nrow(dat) == 2L) {
        "mr_ivw_mre"
      } else {
        c("mr_ivw_mre", "mr_weighted_median", "mr_egger_regression")
      }
      mr_result <- safe_analysis(suppressMessages(mr(dat, method_list = methods)))
      if (mr_result$status == "estimated") {
        estimates <- mr_result$value
        estimates$gwas_accession <- accession
        estimates$taxon_name <- species_name
        estimates$tier <- tier_name
        estimates$sensitivity_model <- model_name
        estimates$analysis_status <- "estimated"
        estimates$error_message <- NA_character_
        mr_rows <- append_row(mr_rows, estimates)
      } else {
        mr_rows <- append_row(
          mr_rows,
          data.frame(
            gwas_accession = accession,
            taxon_name = species_name,
            tier = tier_name,
            sensitivity_model = model_name,
            nsnp = nrow(dat),
            analysis_status = mr_result$status,
            error_message = mr_result$message,
            stringsAsFactors = FALSE
          )
        )
      }

      if (nrow(dat) >= 3L && requireNamespace("mr.raps", quietly = TRUE)) {
        raps_result <- safe_analysis(suppressWarnings(mr.raps::mr.raps(
          b_exp = dat$beta.exposure,
          b_out = dat$beta.outcome,
          se_exp = dat$se.exposure,
          se_out = dat$se.outcome,
          over.dispersion = TRUE,
          loss.function = "huber",
          diagnosis = FALSE
        )))
        if (raps_result$status == "estimated") {
          raps <- raps_result$value
          mr_rows <- append_row(
            mr_rows,
            data.frame(
              id.exposure = unique(dat$id.exposure)[1L],
              id.outcome = unique(dat$id.outcome)[1L],
              outcome = unique(dat$outcome)[1L],
              exposure = unique(dat$exposure)[1L],
              method = "Robust adjusted profile score (RAPS)",
              nsnp = nrow(dat),
              b = raps$beta.hat,
              se = raps$beta.se,
              pval = raps$beta.p.value,
              gwas_accession = accession,
              taxon_name = species_name,
              tier = tier_name,
              sensitivity_model = model_name,
              analysis_status = "estimated",
              error_message = NA_character_,
              stringsAsFactors = FALSE
            )
          )
        } else {
          mr_rows <- append_row(
            mr_rows,
            data.frame(
              gwas_accession = accession,
              taxon_name = species_name,
              tier = tier_name,
              sensitivity_model = model_name,
              nsnp = nrow(dat),
              method = "Robust adjusted profile score (RAPS)",
              analysis_status = raps_result$status,
              error_message = raps_result$message,
              stringsAsFactors = FALSE
            )
          )
        }
      }

      if (nrow(dat) >= 2L) {
        heterogeneity_methods <- "mr_ivw_mre"
        if (nrow(dat) >= 3L) {
          heterogeneity_methods <- c(heterogeneity_methods, "mr_egger_regression")
        }
        heterogeneity_result <- safe_analysis(
          suppressMessages(mr_heterogeneity(dat, method_list = heterogeneity_methods))
        )
        if (heterogeneity_result$status == "estimated") {
          heterogeneity <- heterogeneity_result$value
          heterogeneity$gwas_accession <- accession
          heterogeneity$taxon_name <- species_name
          heterogeneity$tier <- tier_name
          heterogeneity$sensitivity_model <- model_name
          heterogeneity_rows <- append_row(heterogeneity_rows, heterogeneity)
        }
      }

      if (nrow(dat) >= 3L) {
        pleiotropy_result <- safe_analysis(suppressMessages(mr_pleiotropy_test(dat)))
        if (pleiotropy_result$status == "estimated") {
          pleiotropy <- pleiotropy_result$value
          pleiotropy$gwas_accession <- accession
          pleiotropy$taxon_name <- species_name
          pleiotropy$tier <- tier_name
          pleiotropy$sensitivity_model <- model_name
          pleiotropy_rows <- append_row(pleiotropy_rows, pleiotropy)
        }
        leave_one_out_result <- safe_analysis(suppressMessages(mr_leaveoneout(dat)))
        if (leave_one_out_result$status == "estimated") {
          leave_one_out <- leave_one_out_result$value
          leave_one_out$gwas_accession <- accession
          leave_one_out$taxon_name <- species_name
          leave_one_out$tier <- tier_name
          leave_one_out$sensitivity_model <- model_name
          leave_one_out_rows <- append_row(leave_one_out_rows, leave_one_out)
        }
      }

      dat$r.exposure <- get_r_from_lor(
        dat$beta.exposure,
        dat$eaf.exposure,
        dat$ncase.exposure,
        dat$ncontrol.exposure,
        dat$prevalence.exposure
      )
      dat$r.outcome <- get_r_from_lor(
        dat$beta.outcome,
        dat$eaf.outcome,
        dat$ncase.outcome,
        dat$ncontrol.outcome,
        dat$prevalence.outcome
      )
      steiger_result <- safe_analysis(suppressMessages(directionality_test(dat)))
      if (steiger_result$status == "estimated") {
        steiger <- steiger_result$value
        steiger$gwas_accession <- accession
        steiger$taxon_name <- species_name
        steiger$tier <- tier_name
        steiger$sensitivity_model <- model_name
        steiger_rows <- append_row(steiger_rows, steiger)
      }

      if (run_presso && nrow(dat) >= 4L && requireNamespace("MRPRESSO", quietly = TRUE)) {
        presso_seed <- derive_mr_presso_seed(accession, tier_name, model_name)
        set.seed(
          presso_seed,
          kind = "Mersenne-Twister",
          normal.kind = "Inversion",
          sample.kind = "Rejection"
        )
        presso_rng_kind <- paste(RNGkind(), collapse = ";")
        presso_result <- safe_analysis(
          suppressMessages(MRPRESSO::mr_presso(
            BetaOutcome = "beta.outcome",
            BetaExposure = "beta.exposure",
            SdOutcome = "se.outcome",
            SdExposure = "se.exposure",
            OUTLIERtest = TRUE,
            DISTORTIONtest = TRUE,
            data = dat,
            NbDistribution = mr_presso_nb_distribution,
            SignifThreshold = 0.05
          ))
        )
        if (presso_result$status == "estimated") {
          presso <- presso_result$value
          global <- presso[["MR-PRESSO results"]][["Global Test"]]
          outlier <- presso[["MR-PRESSO results"]][["Outlier Test"]]
          distortion <- presso[["MR-PRESSO results"]][["Distortion Test"]]
          presso_rows <- append_row(
            presso_rows,
            data.frame(
              gwas_accession = accession,
              taxon_name = species_name,
              tier = tier_name,
              sensitivity_model = model_name,
              nsnp = nrow(dat),
              rng_seed = presso_seed,
              rng_kind = presso_rng_kind,
              nb_distribution = mr_presso_nb_distribution,
              global_rss_observed = global$RSSobs,
              global_pvalue = as.character(global$Pvalue),
              outlier_count = if (is.null(outlier)) 0L else nrow(outlier),
              distortion_pvalue = if (is.null(distortion)) NA_character_ else as.character(distortion$Pvalue),
              analysis_status = "estimated",
              error_message = NA_character_,
              stringsAsFactors = FALSE
            )
          )
        } else {
          presso_rows <- append_row(
            presso_rows,
            data.frame(
              gwas_accession = accession,
              taxon_name = species_name,
              tier = tier_name,
              sensitivity_model = model_name,
              nsnp = nrow(dat),
              rng_seed = presso_seed,
              rng_kind = presso_rng_kind,
              nb_distribution = mr_presso_nb_distribution,
              analysis_status = presso_result$status,
              error_message = presso_result$message,
              stringsAsFactors = FALSE
            )
          )
        }
      }
    }
  }
}

status_table <- bind_rows_base(analysis_status)
harmonised_table <- bind_rows_base(harmonised_rows)
mr_table <- bind_rows_base(mr_rows)
heterogeneity_table <- bind_rows_base(heterogeneity_rows)
pleiotropy_table <- bind_rows_base(pleiotropy_rows)
steiger_table <- bind_rows_base(steiger_rows)
leave_one_out_table <- bind_rows_base(leave_one_out_rows)
presso_table <- bind_rows_base(presso_rows)

if (nrow(mr_table) > 0L) {
  mr_table$primary_estimator <- with(
    mr_table,
    (nsnp == 1L & method == "Wald ratio") |
      (nsnp >= 2L & method == "Inverse variance weighted (multiplicative random effects)")
  )
  mr_table$nominal_p_lt_0_05 <- !is.na(mr_table$pval) & mr_table$pval < 0.05
  mr_table$odds_ratio <- exp(mr_table$b)
  mr_table$or_ci_lower <- exp(mr_table$b - 1.96 * mr_table$se)
  mr_table$or_ci_upper <- exp(mr_table$b + 1.96 * mr_table$se)
  mr_table$bh_q_primary <- NA_real_
  full_primary <- mr_table$sensitivity_model == "full" & mr_table$primary_estimator %in% TRUE
  for (tier_name in names(tiers)) {
    index <- which(full_primary & mr_table$tier == tier_name & !is.na(mr_table$pval))
    if (length(index) > 0L) {
      mr_table$bh_q_primary[index] <- p.adjust(mr_table$pval[index], method = "BH")
    }
  }
  original_beta <- setNames(crosswalk$mibiogen_primary_beta, crosswalk$accession)
  mr_table$mibiogen_primary_beta <- original_beta[mr_table$gwas_accession]
  mr_table$direction_concordant_with_genus <- sign(mr_table$b) == sign(mr_table$mibiogen_primary_beta)
}

full_primary_table <- mr_table[
  mr_table$sensitivity_model == "full" & mr_table$primary_estimator %in% TRUE,
  ,
  drop = FALSE
]
nominal_candidates <- mr_table[
  mr_table$nominal_p_lt_0_05 %in% TRUE,
  ,
  drop = FALSE
]

write.csv(status_table, file.path(results_dir, "species_mr_analysis_status.csv"), row.names = FALSE)
write.csv(harmonised_table, file.path(results_dir, "species_harmonised_instruments.csv"), row.names = FALSE)
write.csv(mr_table, file.path(results_dir, "species_mr_all_estimators.csv"), row.names = FALSE)
write.csv(full_primary_table, file.path(results_dir, "species_mr_primary_estimates.csv"), row.names = FALSE)
write.csv(nominal_candidates, file.path(results_dir, "species_mr_nominal_p_lt_0_05.csv"), row.names = FALSE)
write.csv(heterogeneity_table, file.path(results_dir, "species_mr_heterogeneity.csv"), row.names = FALSE)
write.csv(pleiotropy_table, file.path(results_dir, "species_mr_egger_intercept.csv"), row.names = FALSE)
write.csv(steiger_table, file.path(results_dir, "species_mr_steiger_directionality.csv"), row.names = FALSE)
write.csv(leave_one_out_table, file.path(results_dir, "species_mr_leave_one_out.csv"), row.names = FALSE)
if (nrow(presso_table) > 0L || !file.exists(presso_output_path)) {
  write.csv(presso_table, presso_output_path, row.names = FALSE)
}

sink(file.path(results_dir, "r_session_info.txt"))
sessionInfo()
sink()

cat("Species-level Mendelian randomization completed:", results_dir, "\n")
