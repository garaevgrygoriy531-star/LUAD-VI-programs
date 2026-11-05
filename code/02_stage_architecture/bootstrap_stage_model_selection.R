# Patient-cluster bootstrap for the gradual, step and hinge stage models.
suppressPackageStartupMessages(library(lme4))
args <- commandArgs(trailingOnly=TRUE)
root <- if(length(args)) args[[1]] else '.'
B <- 1000L
seed <- 20261027L
fit_models <- function(d, levels) {
  d$stage_order <- match(d$pathology, levels)-1L
  span <- length(levels)-2L
  d$pre_order <- pmin(d$stage_order,span)
  d$invasive <- as.integer(d$pathology=='LUAD')
  forms <- list(Gradual=Y~stage_order+(1|patient_id),
                `Invasive step`=Y~invasive+(1|patient_id),
                Hinge=Y~pre_order+invasive+(1|patient_id))
  lapply(forms,function(form) {
    m <- suppressWarnings(lmer(form,data=d,REML=FALSE,
      control=lmerControl(optimizer='bobyqa',optCtrl=list(maxfun=10000),check.conv.singular='ignore')))
    k <- attr(logLik(m),'df'); n <- nobs(m)
    as.numeric(AIC(m)+2*k*(k+1)/(n-k-1))
  })
}
bootstrap_cohort <- function(d, dataset, levels) {
  ids <- unique(as.character(d$patient_id))
  by_patient <- setNames(lapply(ids,function(id)which(as.character(d$patient_id)==id)),ids)
  rows <- vector('list',B)
  for (b in seq_len(B)) {
    set.seed(seed + b + if(dataset=='GSE319666')10000L else 0L)
    draw <- sample(ids,length(ids),replace=TRUE)
    ix <- unlist(by_patient[draw],use.names=FALSE)
    z <- d[ix,,drop=FALSE]
    z$patient_id <- rep(seq_along(draw),lengths(by_patient[draw]))
    scores <- list()
    for (program in c('V2','V3','V4')) {
      z$Y <- z[[program]]
      aicc <- fit_models(z,levels)
      scores[[program]] <- data.frame(replicate=b,program=program,model=names(aicc),AICc=unlist(aicc))
    }
    rows[[b]] <- do.call(rbind,scores)
  }
  do.call(rbind,rows)
}
d1 <- read.delim(file.path(root,'data_processed/pathology/gse308103_stage_scores.tsv'))
d1$pathology <- d1$stage
d2 <- read.delim(file.path(root,'data_processed/pathology/gse319666_sample_scores.tsv'))
d2$pathology <- d2$source_label
d2$pathology[d2$pathology=='N'] <- 'Normal'
d2$pathology[d2$pathology=='ADC'] <- 'LUAD'
all <- rbind(transform(bootstrap_cohort(d1,'GSE308103',c('Normal','AAH','AIS','MIA','LUAD')),dataset='GSE308103'),
             transform(bootstrap_cohort(d2,'GSE319666',c('Normal','AAH','AIS/MIA','LUAD')),dataset='GSE319666'))
all$selected <- ave(all$AICc,all$dataset,all$program,all$replicate,FUN=function(x)as.integer(x==min(x,na.rm=TRUE)))
out <- aggregate(selected~dataset+program+model,all,mean)
dir.create(file.path(root,'results/stage_architecture'),recursive=TRUE,showWarnings=FALSE)
write.csv(all,file.path(root,'results/stage_architecture/model_selection_bootstrap_recomputed_replicates.csv'),row.names=FALSE)
write.csv(out,file.path(root,'results/stage_architecture/model_selection_bootstrap_recomputed_summary.csv'),row.names=FALSE)
