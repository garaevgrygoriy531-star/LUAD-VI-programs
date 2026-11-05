suppressPackageStartupMessages({library(nlme)})
args <- commandArgs(trailingOnly=TRUE)
root <- if(length(args)) args[[1]] else '.'
d <- read.delim(file.path(root,'data_processed/pathology/gse319666_sample_scores.tsv'),check.names=FALSE)
d$in_stageIA_analysis <- tolower(as.character(d$in_stageIA_analysis)) == 'true'
d <- d[d$in_stageIA_analysis,,drop=FALSE]
d$stageIA_group <- factor(d$stageIA_group,levels=c('PML','LUAD_IA'))
stopifnot(nrow(d)==92L,length(unique(d$patient_id))==33L,all(!is.na(d$stageIA_group)))
rows <- lapply(c('V2','V3','V4'),function(program){
  d$score <- d[[program]]
  fit <- lme(score~stageIA_group,random=~1|patient_id,data=d,method='REML')
  tt <- summary(fit)$tTable['stageIA_groupLUAD_IA',]
  critical <- qt(.975,tt['DF'])
  data.frame(program=program,estimate=tt['Value'],ci_low=tt['Value']-critical*tt['Std.Error'],ci_high=tt['Value']+critical*tt['Std.Error'],n_specimens=nrow(d),n_patients=length(unique(d$patient_id)))
})
out <- do.call(rbind,rows)
dir.create(file.path(root,'results/stage_architecture'),recursive=TRUE,showWarnings=FALSE)
write.csv(out,file.path(root,'results/stage_architecture/stage_ia_contrasts.csv'),row.names=FALSE)
print(out)
