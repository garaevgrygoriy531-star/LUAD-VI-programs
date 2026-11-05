#!/usr/bin/env Rscript
# Refit the sample-level pathology contrasts shown for GSE166720 and GSE282617.
# HC3 standard errors are used; confidence intervals use a t reference distribution.
args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args)) args[[1]] else "."
outdir <- file.path(root, "results", "cross_context")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

hc3_contrast <- function(d, score_col, reference, invasive, dataset, program,
                         patient_col = NULL, standardize = FALSE, unit = "sample") {
  d <- d[d$group %in% c(reference, invasive) & is.finite(d[[score_col]]), , drop = FALSE]
  d$group <- factor(d$group, levels = c(reference, invasive))
  if (standardize) d$score <- as.numeric(scale(d[[score_col]])) else d$score <- d[[score_col]]
  fit <- lm(score ~ group, data = d)
  X <- model.matrix(fit)
  residual <- residuals(fit)
  leverage <- hatvalues(fit)
  bread <- solve(crossprod(X))
  weighted_X <- X * (residual / (1 - leverage))^2
  covariance <- bread %*% crossprod(X, weighted_X) %*% bread
  beta <- unname(coef(fit)[2])
  se <- sqrt(covariance[2, 2])
  df <- nrow(d) - ncol(X)
  critical <- qt(.975, df)
  p <- 2 * pt(-abs(beta / se), df)
  n_patients <- if (is.null(patient_col)) NA_integer_ else length(unique(d[[patient_col]]))
  data.frame(dataset = dataset, program = program, reference_group = reference,
             invasive_group = invasive, n_reference = sum(d$group == reference),
             n_invasive = sum(d$group == invasive), n_samples = nrow(d),
             n_patients = n_patients, estimate = beta, se_HC3 = se,
             ci_low = beta - critical * se, ci_high = beta + critical * se,
             p_value = p, inferential_unit = unit,
             model = "Linear contrast with HC3 standard error",
             stringsAsFactors = FALSE)
}

# GSE166720: AIS/MIA versus AC, MP, PAP, and SOL; LPA and unclassified samples are excluded.
g166 <- read.csv(file.path(root, "data_processed", "cross_context",
                           "gse166720_ct_histology_scores.csv"),
                 check.names = FALSE, stringsAsFactors = FALSE)
g166$included_primary <- tolower(as.character(g166$included_primary)) == "true"
g166 <- g166[g166$included_primary, , drop = FALSE]
g166$group <- ifelse(g166$source_label %in% c("AIS", "MIA"), "AIS/MIA", "Invasive")
r166 <- lapply(c("V2", "V3", "V4"), function(p) {
  hc3_contrast(g166, paste0("VI_cluster_", match(p, c("V1", "V2", "V3", "V4")) , "_score"),
               "AIS/MIA", "Invasive", "GSE166720", p, "patient_id",
               standardize = TRUE, unit = "patient")
})

# GSE282617: the public sample manifest does not include a patient crosswalk.
g282 <- read.csv(file.path(root, "data_processed", "cross_context",
                           "gse282617_sample_scores.csv"),
                 check.names = FALSE, stringsAsFactors = FALSE)
g282$in_primary_contrast <- tolower(as.character(g282$in_primary_contrast)) == "true"
g282 <- g282[g282$in_primary_contrast, , drop = FALSE]
g282$group <- ifelse(g282$pathology_group %in% c("AIS", "MIA"), "AIS/MIA", "Invasive")
r282 <- lapply(c("V2", "V3", "V4"), function(p) {
  hc3_contrast(g282, p, "AIS/MIA", "Invasive", "GSE282617", p,
               standardize = FALSE, unit = "sample")
})
out <- do.call(rbind, c(r166, r282))
write.csv(out, file.path(outdir, "pathology_contrasts.csv"), row.names = FALSE, na = "")
message("Cross-context pathology contrasts written to ", outdir)
