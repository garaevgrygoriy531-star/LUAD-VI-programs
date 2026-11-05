suppressPackageStartupMessages({library(data.table);library(lme4);library(lmerTest)})
args <- commandArgs(trailingOnly=TRUE)
root <- if(length(args)) args[[1]] else '.'
outdir <- file.path(root,'results','cellular')
dir.create(outdir,recursive=TRUE,showWarnings=FALSE)
d <- fread(file.path(root,'data_processed','cellular','fine_type_scores.tsv'),data.table=FALSE)
d <- d[d$cell_count>=50 & d$pathology%in%c('AAH','AIS','MIA','LUAD'),,drop=FALSE]
d$LUAD_TRANSITION <- as.integer(d$pathology=='LUAD')
lineages <- c('epithelial','fibroblast_stromal','endothelial','myeloid')
programs <- c('V2','V3','V4')
fit_one <- function(x,program,adjusted=FALSE) {
  x <- x[is.finite(x[[program]]) & !is.na(x$patient_id) & !is.na(x$sample_id),,drop=FALSE]
  precursor <- length(unique(x$patient_id[x$pathology%in%c('AAH','AIS','MIA')]))
  invasive <- length(unique(x$patient_id[x$pathology=='LUAD']))
  if(nrow(x)<8 || precursor<5 || invasive<5) return(NULL)
  x$fine_type <- factor(x$fine_type);x$patient_id <- factor(x$patient_id);x$sample_id <- factor(x$sample_id)
  x$Y <- x[[program]]
  formula <- if(adjusted) Y~LUAD_TRANSITION+V1+fine_type+(1|patient_id)+(1|sample_id) else Y~LUAD_TRANSITION+fine_type+(1|patient_id)+(1|sample_id)
  model <- suppressWarnings(tryCatch(lmerTest::lmer(formula,data=x,REML=FALSE,control=lmerControl(optimizer='bobyqa',optCtrl=list(maxfun=10000),check.conv.singular='ignore')),error=function(e)NULL))
  fallback <- FALSE
  if(is.null(model)||isSingular(model,tol=1e-4)) {
    fallback <- TRUE
    formula <- if(adjusted) Y~LUAD_TRANSITION+V1+fine_type+(1|patient_id) else Y~LUAD_TRANSITION+fine_type+(1|patient_id)
    model <- suppressWarnings(tryCatch(lmerTest::lmer(formula,data=x,REML=FALSE,control=lmerControl(optimizer='bobyqa',optCtrl=list(maxfun=10000),check.conv.singular='ignore')),error=function(e)NULL))
  }
  if(is.null(model)||!'LUAD_TRANSITION'%in%names(fixef(model))) return(NULL)
  sm <- coef(summary(model))['LUAD_TRANSITION',]
  b <- unname(fixef(model)['LUAD_TRANSITION']);df <- sm['df'];se <- sm['Std. Error']
  data.frame(program=program,adjusted_for_V1=adjusted,beta=b,se=se,df=df,ci_low=b-qt(.975,df)*se,ci_high=b+qt(.975,df)*se,p=sm['Pr(>|t|)'],n_units=nrow(x),n_patients=length(unique(x$patient_id)),n_fine_types=nlevels(x$fine_type),patient_only_fallback=fallback)
}
results <- list()
for(lineage in lineages) {
  x <- d[d$broad_lineage==lineage,,drop=FALSE]
  if(!nrow(x)) next
  counts <- lapply(split(x,x$fine_type),function(z)c(precursor=length(unique(z$patient_id[z$pathology%in%c('AAH','AIS','MIA')])),invasive=length(unique(z$patient_id[z$pathology=='LUAD'])),total=length(unique(z$patient_id))))
  eligible <- names(Filter(function(z)z['precursor']>=5&&z['invasive']>=5&&z['total']>=12,counts))
  x <- x[x$fine_type%in%eligible,,drop=FALSE]
  for(program in programs) for(adjusted in c(FALSE,TRUE)) {
    fit <- fit_one(x,program,adjusted)
    if(!is.null(fit)){fit$broad_lineage<-lineage;results[[length(results)+1]]<-fit}
  }
}
if(!length(results)) stop('No eligible fine-cell-type models were estimable from the included processed scores.')
out <- rbindlist(results,fill=TRUE)
out[,q_value:=p.adjust(p,'BH'),by=adjusted_for_V1]
fwrite(out,file.path(outdir,'fine_celltype_models_from_processed_scores.csv'))
print(out)
