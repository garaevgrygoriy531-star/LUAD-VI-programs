#!/usr/bin/env Rscript
# Fit stage-architecture models in two LUAD cohorts.
# Inputs: processed, sample-level stage-score files in data_processed/pathology/.
# Outputs: model estimates, patient-cluster bootstrap intervals, and sensitivity results.
# Inferential unit: patient; specimens are resampled within patient clusters.
# Required package: lme4.

options(stringsAsFactors = FALSE, warn = 1)
suppressPackageStartupMessages(library(lme4))

args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args)) args[[1]] else "."

run_stage_architecture <- function(root) {
  output_dir <- file.path(root, "results", "stage_architecture")
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

  programs <- c("V2", "V3", "V4")
  bootstrap_replicates <- 5000L
  model_forms <- list(
    gradual = Y ~ stage_order + (1 | patient_id),
    step = Y ~ invasive + (1 | patient_id),
    hinge = Y ~ pre_order + invasive + (1 | patient_id)
  )

  fit_model <- function(formula, data) {
    model <- tryCatch(
      suppressWarnings(lmer(
        formula, data = data, REML = FALSE,
        control = lmerControl(
          optimizer = "bobyqa", optCtrl = list(maxfun = 10000),
          check.conv.singular = "ignore"
        )
      )),
      error = function(e) NULL
    )
    if (is.null(model) || any(!is.finite(fixef(model)))) return(NULL)
    model
  }

  calculate_aicc <- function(model) {
    if (is.null(model)) return(NA_real_)
    parameter_count <- attr(logLik(model), "df")
    observations <- nobs(model)
    as.numeric(AIC(model) + 2 * parameter_count * (parameter_count + 1) /
      (observations - parameter_count - 1))
  }

  extract_effects <- function(model, pre_span, direction) {
    missing <- c(
      beta_pre = NA_real_, PRE_SPAN = NA_real_, INV_INCREMENT = NA_real_,
      ALIGNED_PRE_SPAN = NA_real_, ALIGNED_INV_INCREMENT = NA_real_,
      STEP_MINUS_SPAN_CONTRAST = NA_real_, INV_FRACTION = NA_real_
    )
    if (is.null(model)) return(missing)
    coefficients <- fixef(model)
    pre_slope <- unname(coefficients["pre_order"])
    invasion_increment <- unname(coefficients["invasive"])
    pre_change <- pre_span * pre_slope
    c(
      beta_pre = pre_slope,
      PRE_SPAN = pre_change,
      INV_INCREMENT = invasion_increment,
      ALIGNED_PRE_SPAN = direction * pre_change,
      ALIGNED_INV_INCREMENT = direction * invasion_increment,
      STEP_MINUS_SPAN_CONTRAST = direction * (invasion_increment - pre_change),
      INV_FRACTION = abs(invasion_increment) /
        (abs(invasion_increment) + abs(pre_change))
    )
  }

  bootstrap_one <- function(replicate_id) {
    set.seed(seed_base + replicate_id)
    sampled_patients <- sample(patient_ids, length(patient_ids), replace = TRUE)
    row_indices <- unlist(patient_rows[sampled_patients], use.names = FALSE)
    resampled <- data[row_indices, , drop = FALSE]
    resampled$patient_id <- rep(seq_along(sampled_patients),
                                lengths(patient_rows[sampled_patients]))

    do.call(rbind, lapply(programs, function(program) {
      resampled$Y <- resampled[[program]]
      model <- fit_model(model_forms$hinge, resampled)
      effects <- extract_effects(model, pre_span,
                                 if (program == "V4") -1 else 1)
      data.frame(
        replicate = replicate_id, program = program, as.list(effects),
        fit_ok = !is.null(model),
        singular = if (is.null(model)) NA else isSingular(model),
        convergence_message = if (is.null(model)) "FIT_FAILED" else
          paste(model@optinfo$conv$lme4$messages, collapse = ";"),
        stringsAsFactors = FALSE
      )
    }))
  }

  worker_environment <- new.env(parent = .GlobalEnv)
  worker_environment$model_forms <- model_forms
  worker_environment$programs <- programs
  worker_environment$fit_model <- fit_model
  worker_environment$extract_effects <- extract_effects
  worker_environment$bootstrap_one <- bootstrap_one
  environment(worker_environment$fit_model) <- worker_environment
  environment(worker_environment$extract_effects) <- worker_environment
  environment(worker_environment$bootstrap_one) <- worker_environment

  stage_scores_1 <- read.delim(file.path(
    root, "data_processed", "pathology", "gse308103_stage_scores.tsv"
  ))
  stage_scores_2 <- read.delim(file.path(
    root, "data_processed", "pathology", "gse319666_sample_scores.tsv"
  ))
  stage_scores_1$pathology <- stage_scores_1$stage
  stage_scores_2$pathology <- stage_scores_2$source_label
  stage_scores_2$pathology[stage_scores_2$pathology == "N"] <- "Normal"
  stage_scores_2$pathology[stage_scores_2$pathology == "ADC"] <- "LUAD"
  cohorts <- list(GSE308103 = stage_scores_1, GSE319666 = stage_scores_2)

  primary_results <- list()
  patient_exclusion_results <- list()
  lesion_only_results <- list()
  bootstrap_results <- list()
  stage_counts <- list()
  worker_count <- min(10L, max(1L, parallel::detectCores() - 2L))
  cluster <- parallel::makePSOCKcluster(worker_count, rscript_args = "--vanilla")
  on.exit(parallel::stopCluster(cluster), add = TRUE)
  parallel::clusterEvalQ(cluster, {
    options(warn = -1)
    suppressPackageStartupMessages(library(lme4))
    Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1",
               MKL_NUM_THREADS = "1")
    NULL
  })
  for (cohort_index in seq_along(cohorts)) {
    for (scenario_index in 1:2) {
      dataset <- names(cohorts)[cohort_index]
      scenario <- c("PRIMARY_WITH_NORMAL", "LESION_ONLY")[scenario_index]
      data <- cohorts[[cohort_index]]
      if (scenario_index == 2) data <- data[data$pathology != "Normal", , drop = FALSE]
      stage_levels <- if (dataset == "GSE308103") {
        c("Normal", "AAH", "AIS", "MIA", "LUAD")
      } else {
        c("Normal", "AAH", "AIS/MIA", "LUAD")
      }
      if (scenario_index == 2) stage_levels <- stage_levels[stage_levels != "Normal"]
      data$stage_order <- match(data$pathology, stage_levels) - 1
      pre_span <- length(stage_levels) - 2
      data$pre_order <- pmin(data$stage_order, pre_span)
      data$invasive <- as.integer(data$pathology == "LUAD")
      stopifnot(all(is.finite(data$stage_order)))

      patient_ids <- unique(as.character(data$patient_id))
      patient_rows <- setNames(lapply(patient_ids, function(id) {
        which(as.character(data$patient_id) == id)
      }), patient_ids)
      seed_base <- 20261027L + cohort_index * 100000L +
        scenario_index * 10000L
      worker_environment$data <- data
      worker_environment$pre_span <- pre_span
      worker_environment$patient_ids <- patient_ids
      worker_environment$patient_rows <- patient_rows
      worker_environment$seed_base <- seed_base
      observed <- list()
      for (program in programs) {
        model_data <- data
        model_data$Y <- model_data[[program]]
        hinge_model <- fit_model(model_forms$hinge, model_data)
        stopifnot(!is.null(hinge_model))
        effects <- extract_effects(hinge_model, pre_span,
                                   if (program == "V4") -1 else 1)
        observed[[program]] <- list(
          effects = effects, model = hinge_model,
          AICc_hinge = calculate_aicc(hinge_model),
          AICc_step = calculate_aicc(fit_model(model_forms$step, model_data)),
          AICc_gradual = calculate_aicc(fit_model(model_forms$gradual, model_data))
        )

        for (patient in patient_ids) {
          leave_one_out <- model_data[
            as.character(model_data$patient_id) != patient, , drop = FALSE
          ]
          leave_one_model <- fit_model(model_forms$hinge, leave_one_out)
          leave_one_effects <- extract_effects(leave_one_model, pre_span,
                                               if (program == "V4") -1 else 1)
          patient_exclusion_results[[length(patient_exclusion_results) + 1]] <-
            data.frame(
              dataset = dataset, scenario = scenario, program = program,
              removed_patient = patient, as.list(leave_one_effects),
              step_minus_span_sign_changed = if (is.null(leave_one_model)) NA else
                sign(leave_one_effects["STEP_MINUS_SPAN_CONTRAST"]) !=
                  sign(effects["STEP_MINUS_SPAN_CONTRAST"]),
              fit_ok = !is.null(leave_one_model), stringsAsFactors = FALSE
            )
        }
      }

      message("Running patient-cluster bootstrap: ", dataset, " / ", scenario,
              " (", length(patient_ids), " patients; ", nrow(data), " specimens).")
      replicate_chunks <- split(seq_len(bootstrap_replicates),
                                ceiling(seq_len(bootstrap_replicates) / 250L))
      chunk_results <- vector("list", length(replicate_chunks))
      for (chunk_index in seq_along(replicate_chunks)) {
        chunk_results[[chunk_index]] <- do.call(rbind, parallel::parLapply(
          cluster, replicate_chunks[[chunk_index]], worker_environment$bootstrap_one
        ))
      }
      bootstrap <- do.call(rbind, chunk_results)
      bootstrap$dataset <- dataset
      bootstrap$scenario <- scenario
      bootstrap$seed_base <- seed_base
      write.csv(bootstrap, file.path(output_dir,
        paste0("BOOTSTRAP_", dataset, "_", scenario, ".csv")), row.names = FALSE)
      bootstrap_results[[paste(dataset, scenario)]] <- bootstrap

      for (program in programs) {
        model_result <- observed[[program]]
        effects <- model_result$effects
        program_bootstrap <- bootstrap[bootstrap$program == program, ]
        result <- data.frame(
          dataset = dataset, scenario = scenario, program = program,
          n_specimens = nrow(data), n_patients = length(patient_ids),
          pre_span_categories = pre_span, as.list(effects),
          AICc_hinge = model_result$AICc_hinge,
          AICc_step = model_result$AICc_step,
          AICc_gradual = model_result$AICc_gradual,
          bootstrap_requested = bootstrap_replicates,
          bootstrap_valid = sum(program_bootstrap$fit_ok),
          bootstrap_singular = sum(program_bootstrap$singular, na.rm = TRUE),
          seed_base = seed_base,
          estimator = "Linear mixed-effects model (maximum likelihood; bobyqa optimizer)",
          stringsAsFactors = FALSE
        )
        for (effect_name in names(effects)) {
          interval <- quantile(
            program_bootstrap[[effect_name]][is.finite(program_bootstrap[[effect_name]])],
            c(0.025, 0.975), names = FALSE
          )
          result[[paste0(effect_name, "_ci_low")]] <- interval[1]
          result[[paste0(effect_name, "_ci_high")]] <- interval[2]
        }
        if (scenario_index == 1) {
          primary_results[[length(primary_results) + 1]] <- result
        } else {
          lesion_only_results[[length(lesion_only_results) + 1]] <- result
        }
      }

      for (stage in stage_levels) {
        stage_order <- match(stage, stage_levels) - 1
        stage_counts[[length(stage_counts) + 1]] <- data.frame(
          dataset = dataset, scenario = scenario, pathology = stage,
          stage_order = stage_order, pre_order = min(stage_order, pre_span),
          LUAD_indicator = as.integer(stage == "LUAD"),
          n_specimens = sum(data$pathology == stage),
          n_patients = length(unique(data$patient_id[data$pathology == stage])),
          stringsAsFactors = FALSE
        )
      }
      write.csv(do.call(rbind, primary_results),
                file.path(output_dir, "primary_estimates.csv"), row.names = FALSE)
      if (length(lesion_only_results)) {
        write.csv(do.call(rbind, lesion_only_results),
                  file.path(output_dir, "lesion_only_sensitivity.csv"), row.names = FALSE)
      }
      write.csv(do.call(rbind, patient_exclusion_results),
                file.path(output_dir, "leave_one_patient_out.csv"), row.names = FALSE)
    }
  }

  write.csv(do.call(rbind, stage_counts), file.path(output_dir, "stage_coding.csv"),
            row.names = FALSE)
  saveRDS(do.call(rbind, bootstrap_results),
          file.path(output_dir, "bootstrap_distributions.rds"))
  message("Stage architecture analysis complete.")
}

run_stage_architecture(root)
