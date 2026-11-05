# Refit the Figure 3A contrasts with patient-cluster bootstrap confidence intervals.
# Run from any directory with the repository root as the first argument.
args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args) >= 1L) args[[1L]] else "."
out_dir <- file.path(root, "results", "pathology")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
out_file <- file.path(out_dir, "v1_adjusted_stage_models.csv")

suppressPackageStartupMessages(library(lme4))
reps <- 1000L
seed <- 20261016L
outcomes <- c("V2", "V3", "V4")
complete_score_columns <- c("ORACLE", "V2", "V3", "V4", "V1", "E", "patient")
control <- lmerControl(
  optimizer = "bobyqa",
  optCtrl = list(maxfun = 100000),
  check.conv.singular = "ignore"
)

read_cohort <- function(dataset, estimand = "LUAD_transition") {
  if (dataset == "GSE308103") {
    d <- read.delim(
      file.path(root, "data_processed", "pathology", "gse308103_stage_scores.tsv"),
      check.names = FALSE
    )
    d$E <- if (estimand == "any_invasive") {
      as.integer(d$stage %in% c("MIA", "LUAD"))
    } else {
      as.integer(d$stage == "LUAD")
    }
  } else {
    d <- read.delim(
      file.path(root, "data_processed", "pathology", "gse319666_sample_scores.tsv"),
      check.names = FALSE
    )
    d$E <- as.integer(d$source_label == "ADC")
  }
  d$patient <- factor(d$patient_id)
  complete <- complete.cases(d[, complete_score_columns, drop = FALSE])
  d <- d[complete, , drop = FALSE]
  d$patient <- droplevels(d$patient)
  d
}

make_patient_draws <- function(d) {
  patients <- levels(d$patient)
  matrix(sample(patients, reps * length(patients), replace = TRUE),
         nrow = reps, ncol = length(patients))
}

fit_model <- function(formula, data) {
  lmer(formula, data = data, REML = FALSE, control = control)
}

bootstrap_effects <- function(d, draws, outcome, adjusted) {
  formula <- if (adjusted) {
    as.formula(paste(outcome, "~ E + V1 + (1|bootpatient)"))
  } else {
    as.formula(paste(outcome, "~ E + (1|bootpatient)"))
  }
  values <- rep(NA_real_, reps)
  patient_rows <- split(d, d$patient)
  for (i in seq_len(reps)) {
    parts <- vector("list", ncol(draws))
    for (j in seq_len(ncol(draws))) {
      z <- patient_rows[[draws[i, j]]]
      z$bootpatient <- paste0("b", j)
      parts[[j]] <- z
    }
    boot_data <- do.call(rbind, parts)
    values[i] <- tryCatch(
      unname(fixef(fit_model(formula, boot_data))["E"]),
      error = function(e) NA_real_
    )
  }
  values
}

summarize_cohort <- function(dataset, d, draws) {
  rows <- list()
  k <- 0L
  for (outcome in outcomes) {
    for (adjusted in c(FALSE, TRUE)) {
      full_formula <- if (adjusted) {
        as.formula(paste(outcome, "~ E + V1 + (1|patient)"))
      } else {
        as.formula(paste(outcome, "~ E + (1|patient)"))
      }
      boot <- bootstrap_effects(d, draws, outcome, adjusted)
      estimate <- unname(fixef(fit_model(full_formula, d))["E"])
      finite <- boot[is.finite(boot)]
      interval <- if (length(finite)) {
        as.numeric(quantile(finite, probs = c(0.025, 0.975), names = FALSE))
      } else {
        c(NA_real_, NA_real_)
      }
      k <- k + 1L
      rows[[k]] <- data.frame(
        dataset = dataset,
        outcome = outcome,
        adjustment = if (adjusted) "V1-adjusted" else "Unadjusted",
        estimate = estimate,
        bootstrap_ci_low = interval[1],
        bootstrap_ci_high = interval[2],
        model_formula = paste(deparse(full_formula), collapse = ""),
        n_specimens = nrow(d),
        n_patients = nlevels(d$patient),
        bootstrap_replicates = reps,
        bootstrap_success = length(finite),
        seed = seed,
        optimizer = "bobyqa",
        estimation = "maximum likelihood",
        stringsAsFactors = FALSE
      )
    }
  }
  do.call(rbind, rows)
}

set.seed(seed)
# Draw order follows the two pathology contrasts in GSE308103 and then GSE319666.
gse308103 <- read_cohort("GSE308103", "LUAD_transition")
draws_308103 <- make_patient_draws(gse308103)
gse308103_any_invasive <- read_cohort("GSE308103", "any_invasive")
invisible(make_patient_draws(gse308103_any_invasive))
gse319666 <- read_cohort("GSE319666")
draws_319666 <- make_patient_draws(gse319666)

result <- rbind(
  summarize_cohort("GSE308103", gse308103, draws_308103),
  summarize_cohort("GSE319666", gse319666, draws_319666)
)
write.csv(result, out_file, row.names = FALSE, na = "")
cat("Wrote", nrow(result), "model summaries to", out_file,
    "with", reps, "patient-cluster bootstrap replicates per contrast.\n")