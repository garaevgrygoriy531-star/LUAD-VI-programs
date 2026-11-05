#!/usr/bin/env Rscript
# Refit primary stage I overall-survival associations from the included cohort table.
# Statistical unit: patient. Scores are standardized within each cohort.
suppressPackageStartupMessages({
  library(data.table)
  library(survival)
})
args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args)) args[[1]] else "."
input <- file.path(root, "data_processed", "outcome", "stageI_primary_cohort_scores.csv")
outdir <- file.path(root, "results", "stage_i_survival")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
d <- fread(input, data.table = FALSE)
cohorts <- c("GSE50081", "GSE31210", "GSE30219")
programs <- c("V2", "V3", "V4")
stopifnot(all(c("cohort", "patient_id", "os_time", "os_event", "primary_eligible", programs) %in% names(d)))
d$primary_eligible <- tolower(as.character(d$primary_eligible)) == "true"
fits <- list()
for (score in programs) for (cohort in cohorts) {
  x <- d[d$primary_eligible & d$cohort == cohort, , drop = FALSE]
  x <- x[is.finite(x$os_time) & x$os_time > 0 & x$os_event %in% c(0, 1) & is.finite(x[[score]]), , drop = FALSE]
  x$score_z <- as.numeric(scale(x[[score]]))
  fit <- coxph(Surv(os_time, os_event) ~ score_z, data = x, ties = "efron")
  z <- summary(fit)$coefficients["score_z", ]
  beta <- unname(z["coef"]); se <- unname(z["se(coef)"])
  fits[[length(fits) + 1L]] <- data.frame(
    cohort = cohort, endpoint = "OS", score = score, model = "PRIMARY_UNADJUSTED",
    n = nrow(x), patients = length(unique(x$patient_id)), events = sum(x$os_event == 1),
    beta = beta, se = se, hr = exp(beta), ci_low = exp(beta - 1.96 * se),
    ci_high = exp(beta + 1.96 * se), p = unname(z["Pr(>|z|)"]),
    formula = "Surv(os_time, os_event) ~ score_z", stringsAsFactors = FALSE
  )
}
primary <- rbindlist(fits)
fwrite(primary, file.path(outdir, "primary_os_cox.csv"))
# Random-effects pooling by restricted maximum likelihood with Hartung-Knapp intervals.
meta_one <- function(y, se) {
  v <- se^2; k <- length(y)
  objective <- function(tau2) {
    w <- 1 / (v + tau2); mu <- sum(w * y) / sum(w)
    sum(log(v + tau2)) + log(sum(w)) + sum(w * (y - mu)^2)
  }
  upper <- max(1, 100 * var(y), 100 * max(v))
  opt <- optimize(objective, c(0, upper))
  tau2 <- if (objective(0) <= opt$objective) 0 else opt$minimum
  w <- 1 / (v + tau2); mu <- sum(w * y) / sum(w)
  q <- sum(w * (y - mu)^2)
  hk_se <- sqrt((q / (k - 1)) / sum(w))
  crit <- qt(.975, df = k - 1)
  data.frame(k = k, beta = mu, se_hk = hk_se, hr = exp(mu),
             ci_low = exp(mu - crit * hk_se), ci_high = exp(mu + crit * hk_se),
             p_hk = 2 * pt(-abs(mu / hk_se), df = k - 1),
             tau2 = tau2, I2 = if (q > 0) max(0, (q - (k - 1)) / q) * 100 else 0)
}
pooled <- rbindlist(lapply(programs, function(s) {
  x <- primary[score == s]
  cbind(data.frame(endpoint = "OS", model = "PRIMARY_UNADJUSTED", score = s,
                   cohorts = paste(x$cohort, collapse = ";")), meta_one(x$beta, x$se))
}))
fwrite(pooled, file.path(outdir, "pooled_os_reml_hk.csv"))
message("Primary stage I OS models written to ", outdir)

