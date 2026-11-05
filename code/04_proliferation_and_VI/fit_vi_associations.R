#!/usr/bin/env Rscript
# Fit prespecified microscopic vascular-invasion models from lesion-level scores.
# Statistical unit: lesion; the duplicate-patient sensitivity retains one lesion at a time.
suppressPackageStartupMessages({ library(logistf); library(brglm2) })
args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args)) args[[1]] else "."
outdir <- file.path(root, "results", "vascular_invasion")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
d <- read.csv(file.path(root, "data_processed", "pathology", "gse273528_vi_samples.csv"),
              check.names = FALSE, stringsAsFactors = FALSE)
stopifnot(nrow(d) == 59L, length(unique(d$patient_id)) == 58L, sum(d$VI == 1) == 17L)
models <- list(VI_V2_plus_V1 = c("V2", "V1"), VI_V3_plus_V1 = c("V3", "V1"), VI_ORACLE = "ORACLE")
fit_set <- function(x, label) {
  for (s in c("V1", "V2", "V3", "ORACLE")) x[[paste0("z", s)]] <- as.numeric(scale(x[[s]]))
  rows <- list()
  for (nm in names(models)) {
    predictor <- paste0("z", models[[nm]][1])
    rhs <- paste0("z", models[[nm]], collapse = " + ")
    form <- as.formula(paste("VI ~", rhs))
    f <- logistf(form, data = x, pl = TRUE,
                 control = logistf.control(maxit = 1000, maxstep = 0.5, lconv = 1e-7))
    j <- match(predictor, names(f$coefficients))
    rows[[length(rows) + 1L]] <- data.frame(
      subset = label, model = nm, predictor = predictor, n = nrow(x),
      patients = length(unique(x$patient_id)), events = sum(x$VI == 1),
      beta = unname(f$coefficients[j]), OR = exp(unname(f$coefficients[j])),
      CI_low = exp(unname(f$ci.lower[j])), CI_high = exp(unname(f$ci.upper[j])),
      p = unname(f$prob[j]), iterations = paste(f$iter, collapse = ";"),
      package = "logistf_profile_PL", stringsAsFactors = FALSE
    )
  }
  do.call(rbind, rows)
}
primary <- fit_set(d, "FULL_59_LESIONS")
primary$BH_q_primary <- NA_real_
ix <- primary$model %in% c("VI_V2_plus_V1", "VI_V3_plus_V1")
primary$BH_q_primary[ix] <- p.adjust(primary$p[ix], method = "BH")
write.csv(primary, file.path(outdir, "primary_firth_models.csv"), row.names = FALSE, na = "NA")
duplicates <- unique(d$patient_id[duplicated(d$patient_id) | duplicated(d$patient_id, fromLast = TRUE)])
sens <- lapply(seq_len(sum(d$patient_id == duplicates[1])), function(i) {
  keep <- d[d$patient_id != duplicates[1] | d$sample_id == d$sample_id[d$patient_id == duplicates[1]][i], ]
  fit_set(keep, paste0("KEEP_", d$sample_id[d$patient_id == duplicates[1]][i]))
})
sensitivity <- do.call(rbind, sens)
write.csv(sensitivity, file.path(outdir, "duplicate_patient_sensitivity.csv"), row.names = FALSE, na = "NA")
message("Microscopic vascular-invasion models written to ", outdir)

