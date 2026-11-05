suppressPackageStartupMessages(library(data.table))
args <- commandArgs(trailingOnly=TRUE)
root <- if(length(args)) args[[1]] else '.'
seed <- 20261027L
B <- 5000L
expected <- c(V1=NA_character_,V2='positive',V3='positive',V4='negative')
exact_permutation_p <- function(y,group) {
  i <- which(group=='IAC');e <- which(group=='AIS_MIA')
  if(!length(i)||!length(e)) return(NA_real_)
  obs <- mean(y[i])-mean(y[e]);combinations <- combn(seq_along(y),length(i))
  differences <- apply(combinations,2,function(ii)mean(y[ii])-mean(y[-ii]))
  mean(abs(differences)>=abs(obs)-1e-12)
}
hedges_g <- function(y,group) {
  i <- y[group=='IAC'];e <- y[group=='AIS_MIA'];ni<-length(i);ne<-length(e)
  if(ni<2L||ne<2L)return(NA_real_)
  pooled<-sqrt(((ni-1)*var(i)+(ne-1)*var(e))/(ni+ne-2))
  if(!is.finite(pooled)||pooled<=0)return(NA_real_)
  (1-3/(4*(ni+ne)-9))*(mean(i)-mean(e))/pooled
}
bootstrap_effects <- function(y,group,seed_local) {
  set.seed(seed_local);i<-y[group=='IAC'];e<-y[group=='AIS_MIA']
  bd<-numeric(B);bg<-numeric(B)
  for(b in seq_len(B)) {ib<-sample(i,length(i),replace=TRUE);eb<-sample(e,length(e),replace=TRUE);bd[b]<-mean(ib)-mean(eb);bg[b]<-hedges_g(c(ib,eb),c(rep('IAC',length(ib)),rep('AIS_MIA',length(eb))))}
  q <- function(x,p) { x<-x[is.finite(x)]; if(!length(x)) NA_real_ else unname(quantile(x,p)) }
  c(diff_low=q(bd,.025),diff_high=q(bd,.975),g_low=q(bg,.025),g_high=q(bg,.975))
}
score_file<-file.path(root,'data_processed','cellular','scrna_pseudobulk_scores.csv')
sc<-fread(score_file,data.table=FALSE)
lineages<-unique(sc$lineage);programs<-c('V1','V2','V3','V4')
primary<-list();influence<-list();er<-0L
for(lineage in lineages) for(program in programs) {
  d<-sc[sc$lineage==lineage&sc$anchor==program&sc$coverage>=.80&is.finite(sc$score)&sc$cell_count>=50,,drop=FALSE]
  d<-d[!duplicated(d$patient_id),,drop=FALSE]
  if(!nrow(d))next
  difference<-mean(d$score[d$group=='IAC'])-mean(d$score[d$group=='AIS_MIA'])
  boot<-bootstrap_effects(d$score,d$group,seed+er+100L)
  direction<-unname(expected[program]);match<-if(is.na(direction))NA else if(direction=='positive')difference>0 else difference<0
  loo<-vapply(d$patient_id,function(pid){z<-d[d$patient_id!=pid,,drop=FALSE];mean(z$score[z$group=='IAC'])-mean(z$score[z$group=='AIS_MIA'])},numeric(1))
  v1<-sc[sc$lineage==lineage&sc$anchor=='V1'&sc$coverage>=.80&is.finite(sc$score)&sc$cell_count>=50,,drop=FALSE];v1<-v1[!duplicated(v1$patient_id),,drop=FALSE]
  paired<-merge(d[,c('patient_id','score','group')],v1[,c('patient_id','score')],by='patient_id',suffixes=c('_target','_V1'))
  v1diff<-NA_real_;v1p<-NA_real_;v1ok<-FALSE
  if(program!='V1'&&nrow(paired)>=5&&length(unique(paired$score_V1))>=3){m<-lm(score_target~score_V1,data=paired);paired$resid<-residuals(m);v1diff<-mean(paired$resid[paired$group=='IAC'])-mean(paired$resid[paired$group=='AIS_MIA']);v1p<-exact_permutation_p(paired$resid,paired$group);v1ok<-TRUE}
  if(program=='V1'){v1diff<-NA_real_;v1p<-NA_real_;v1ok<-NA}
  er<-er+1L
  primary[[er]]<-data.frame(lineage=lineage,anchor=program,n_AIS_MIA=sum(d$group=='AIS_MIA'),n_IAC=sum(d$group=='IAC'),mean_AIS_MIA=mean(d$score[d$group=='AIS_MIA']),mean_IAC=mean(d$score[d$group=='IAC']),difference_IAC_minus_AIS_MIA=difference,hedges_g=hedges_g(d$score,d$group),diff_boot_low=boot['diff_low'],diff_boot_high=boot['diff_high'],hedges_g_boot_low=boot['g_low'],hedges_g_boot_high=boot['g_high'],exact_permutation_p=exact_permutation_p(d$score,d$group),expected_direction=direction,expected_direction_match=match,n_patients=nrow(d),lopo_expected_direction_stability=if(is.na(direction))NA else mean(if(direction=='positive')loo>0 else loo<0),lopo_difference_min=min(loo),lopo_difference_max=max(loo),coverage_min=min(d$coverage),min_cells_per_patient=min(d$cell_count),v1_sensitivity_estimable=v1ok,v1_sensitivity_difference=v1diff,v1_sensitivity_exact_p=v1p,analysis_unit='patient-level lineage pseudobulk; exact stage-label permutation; bootstrap CI')
  influence[[er]]<-data.frame(lineage=lineage,anchor=program,dropped_patient=d$patient_id,full_difference=difference,leave_one_out_difference=loo)
}
sc_effects<-rbindlist(primary,fill=TRUE)
sc_effects$p_adj_BH_by_anchor<-NA_real_
for(program in c('V2','V3','V4')) {ix<-which(sc_effects$anchor==program&is.finite(sc_effects$exact_permutation_p));sc_effects$p_adj_BH_by_anchor[ix]<-p.adjust(sc_effects$exact_permutation_p[ix],method='BH')}
sp<-fread(file.path(root,'data_processed','cellular','spatial_sample_program_scores.csv'),data.table=FALSE)
spout<-list()
for(j in seq_along(programs)) {
  program<-programs[j];d<-sp[sp$anchor==program&sp$coverage>=.80&is.finite(sp$score),,drop=FALSE]
  diff<-mean(d$score[d$group=='IAC'])-mean(d$score[d$group=='AIS_MIA'])
  boot<-bootstrap_effects(d$score,d$group,seed+j+999L)
  direction<-unname(expected[program]);match<-if(is.na(direction))NA else if(direction=='positive')diff>0 else diff<0
  spout[[j]]<-data.frame(anchor=program,n_AIS_MIA=sum(d$group=='AIS_MIA'),n_IAC=sum(d$group=='IAC'),mean_AIS_MIA=mean(d$score[d$group=='AIS_MIA']),mean_IAC=mean(d$score[d$group=='IAC']),difference_IAC_minus_AIS_MIA=diff,hedges_g=hedges_g(d$score,d$group),diff_boot_low=boot['diff_low'],diff_boot_high=boot['diff_high'],hedges_g_boot_low=boot['g_low'],hedges_g_boot_high=boot['g_high'],exact_permutation_p=exact_permutation_p(d$score,d$group),expected_direction=direction,expected_direction_match=match,coverage_min=min(d$coverage),analysis_unit='spatial sample; spots aggregated within sample before normalization')
}
outdir<-file.path(root,'results','cellular');dir.create(outdir,recursive=TRUE,showWarnings=FALSE)
fwrite(sc_effects,file.path(outdir,'gse189357_lineage_effects.csv'))
fwrite(rbindlist(influence,fill=TRUE),file.path(outdir,'gse189357_patient_influence.csv'))
fwrite(rbindlist(spout,fill=TRUE),file.path(outdir,'gse189487_spatial_effects.csv'))
