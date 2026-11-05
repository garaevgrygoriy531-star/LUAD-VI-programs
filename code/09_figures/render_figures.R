suppressPackageStartupMessages({library(ggplot2);library(patchwork);library(grid);library(data.table)})
options(warn=1)
args <- commandArgs(trailingOnly=TRUE)
root <- if(length(args)) args[[1]] else '.'
phase<-if(length(args)>1)args[[2]]else'.'
only_fig<-if(length(args)>2)args[-c(1,2)]else character()
publication<-TRUE
connector_mode<-strsplit(Sys.getenv('LUAD_CONNECTOR_MODE','F3a'),',',fixed=TRUE)[[1]]
datadir <- file.path(root,'data_processed','plot_ready')
outdir<-file.path(root,phase)
dir.create(outdir,recursive=TRUE,showWarnings=FALSE)
for(d in c('PANELS','FIGURES','QC','ACTUAL_SIZE_PROOFS'))dir.create(file.path(outdir,d),recursive=TRUE,showWarnings=FALSE)
# Build grid metrics on the same editable Cairo/Arial device used for export.
metric_file <- tempfile('luad_metric_',fileext='.pdf')
grDevices::cairo_pdf(metric_file,width=10,height=10,family='Arial')
spec <- as.data.frame(fread(file.path(root,'metadata','panel_manifest.csv')))
width_mm = 183
pal <- c(V1='#7D8790',V2='#087F8C',V3='#B17C16',V4='#79559A')
component_pal <- c('Within-cell-type'='#58666F','Fine-cell-type composition'='#D3D8DC','Covariance'='#A69D93','Total'='#222222')
lineage_order <- c('Epithelial','Fibroblast/stromal','Endothelial','Myeloid')
helper_colour<-if(publication)'#E0E3EB'else'#BBBBBB'
component_ink<-c('Within-cell-type'='#58666F','Fine-cell-type composition'='#8C959D','Covariance'='#958A7E','Total'='#222222')
theme_set(theme_classic(base_size=6.5,base_family='Arial')+
  theme(text=element_text(colour='#161616'),axis.text=element_text(colour='#161616',size=6.1),
    axis.title=element_text(size=6.5,face='plain'),axis.line=element_line(linewidth=.22,colour='#51545C'),axis.ticks=element_line(linewidth=.18,colour='#51545C'),
    axis.ticks.length=unit(1,'mm'),legend.title=element_text(size=6),legend.text=element_text(size=5.8),
    legend.key.size=unit(2.5,'mm'),legend.spacing.x=unit(1,'mm'),legend.margin=margin(0,0,0,0),
    strip.background=element_blank(),strip.text=element_text(size=6.4,face='plain'),panel.spacing=unit(1.6,'mm'),
    plot.margin=margin(1,2,1,1,'mm'),plot.title=element_text(size=6.5,face='plain'),plot.caption=element_text(size=5.6,hjust=0),
    legend.position='top'))

rd <- function(id,suffix='') {
  f <- paste0(id,if(nzchar(suffix)) paste0('_',suffix) else '', '.csv')
  d <- as.data.frame(fread(file.path(datadir,f),na.strings=c('','NA','NaN')))
  for(k in intersect(names(d),c('estimate','attenuation_percent','ci_low','ci_high','value','delta','first','second','coverage','pre','luad','match','nes','fdr','number_measured','number_total','coverage_fraction','loo_sign_consistent','loo_significant','loo_total','primary','weight','preference','n','events','patients','samples','gate','mapped','total'))) {
    if(k=='gate') next
    d[[k]] <- suppressWarnings(as.numeric(d[[k]]))
  }
  if('row'%in%names(d))d$row<-gsub('*','',d$row,fixed=TRUE)
  if(id=='S2a' && 'row'%in%names(d)){d$row<-gsub('INV','Invasive',d$row,fixed=TRUE);d<-d[is.finite(d$estimate),,drop=FALSE]}
  if('program'%in%names(d)) d$program <- factor(d$program,levels=c('V1','V2','V3','V4','ORACLE'))
  d
}
order_row <- function(d,levels=NULL) {
  if(is.null(levels)) levels <- if(all(unique(d$row)%in%lineage_order)) lineage_order else unique(d$row)
  d$row <- factor(d$row,levels=rev(levels));d
}
null_line <- function(log=FALSE) geom_vline(xintercept=if(log) 1 else 0,linetype='dashed',linewidth=.2,colour=helper_colour)
ci_geom <- function(mapping=NULL,width=.10,...) geom_errorbar(mapping=mapping,...,orientation='y',width=width,linewidth=.28,na.rm=TRUE,lineend='butt')
base_intervals <- function(d,xlab='Contrast',facets=FALSE,log=FALSE,counts=FALSE) {
  d <- order_row(d)
  if(counts && 'n'%in%names(d)) {
    labels <- as.character(d$row)
    labmap <- setNames(paste0(labels,'  ',d$n,'/',d$events),labels)
    d$row <- factor(as.character(d$row),levels=levels(d$row),labels=labmap[levels(d$row)])
  }
  p <- ggplot(d,aes(x=estimate,y=row,colour=program))+null_line(log)+
    ci_geom(aes(xmin=ci_low,xmax=ci_high))+geom_point(size=1.4,na.rm=TRUE)+
    scale_colour_manual(values=pal,drop=FALSE,guide='none')+labs(x=xlab,y=NULL)
  if(facets)p<-p+facet_wrap(~program,nrow=1)
  if(log) {
    stopifnot(all(d$ci_low[is.finite(d$ci_low)]>0))
    p<-p+scale_x_log10(breaks=c(.1,1,10,100),labels=c('0.1','1','10','100'))
  }
  p
}
models_plot <- function(id,v) {
  d<-rd(id); d$model<-factor(d$model,levels=c('Gradual','Invasive step','Hinge'));d$program<-factor(d$program,levels=c('V4','V3','V2'))
  if(v=='B') {
    d$lab<-sprintf('%.2f',d$delta)
    p<-ggplot(d,aes(model,program,fill=weight))+geom_tile(colour='white',linewidth=.6)+
      geom_text(aes(label=lab),size=2.1,colour='#161616')+scale_fill_gradient(low='#F5F5F2',high='#52798B',limits=c(0,1),breaks=c(0,.5,1),labels=c('0','0.5','1'),name='Akaike weight',guide=guide_colourbar(title.position='top',barwidth=unit(23,'mm'),barheight=unit(2,'mm')))+
      labs(x='Model (cell labels: delta AICc)',y=NULL)+theme(axis.line=element_blank(),axis.ticks=element_blank(),axis.text.x=element_text(angle=25,hjust=1))
    pref<-d[d$model=='Invasive step',];pref$label<-sprintf('Step selection: %.1f%%',pref$preference*100)
    p<-p+geom_text(data=pref,aes(x=3.55,y=program,label=label),inherit.aes=FALSE,size=1.85,hjust=0)+coord_cartesian(xlim=c(.5,5.3),clip='off')
  } else {
    d$y<-as.numeric(d$program)+(as.numeric(d$model)-2)*.18
    p<-ggplot(d,aes(delta,y,colour=program,shape=model))+geom_vline(xintercept=0,linetype='dashed',colour=helper_colour,linewidth=.2)+
      geom_point(size=1.5)+
      scale_colour_manual(values=pal,guide='none')+scale_shape_manual(values=c(1,15,17),name=NULL)+
      scale_y_continuous(breaks=1:3,labels=c('V4','V3','V2'),limits=c(.55,3.6))+labs(x='Delta AICc',y=NULL)+theme(legend.position='bottom',legend.text=element_text(size=5.4),legend.key.width=unit(2.5,'mm'))
    pref<-d[d$model=='Invasive step',]; x<-max(d$delta)+3
    p<-p+geom_text(data=pref,aes(x=x,y=as.numeric(program),label=sprintf('%.1f%%',preference*100)),inherit.aes=FALSE,size=1.85,hjust=0)+
      annotate('text',x=x,y=3.5,label='Step selection',size=1.855,hjust=0)+
      expand_limits(x=x+max(d$delta)*.35)
  }
  p
}
summary_connector <- function(id,d,series_col,open_label,filled_label) {
 if(!id%in%connector_mode)return(NULL)
 a<-d[d[[series_col]]==open_label,];b<-d[d[[series_col]]==filled_label,]
 j<-merge(a,b,by=c('dataset','program'),suffixes=c('_open','_filled'))
 stopifnot(nrow(j)==6)
 geom_segment(data=j,aes(x=estimate_open,xend=estimate_filled,y=y_open,yend=y_filled),inherit.aes=FALSE,colour='#B9BDC2',linewidth=.12,lineend='butt')
}
component_plot <- function(id,v) {
  d<-rd(id,'components');d$component<-factor(d$component,levels=c('Earlier-state span','Invasive-state increment'))
  if(v=='B') {
    ggplot(d,aes(component,estimate,colour=program,group=program))+geom_hline(yintercept=0,colour=helper_colour,linewidth=.2,linetype='dashed')+
      geom_line(linewidth=.35)+geom_errorbar(aes(ymin=ci_low,ymax=ci_high),width=.12,linewidth=.3)+geom_point(aes(shape=component),size=1.5)+
      scale_colour_manual(values=pal,guide='none')+scale_shape_manual(values=c(1,15),guide='none')+facet_grid(dataset~program)+labs(x=NULL,y='Aggressive-direction aligned estimate')+
      theme(axis.text.x=element_text(angle=20,hjust=1,size=5.7))
  } else {
    d$y<-match(as.character(d$program),c('V4','V3','V2'))+ifelse(d$component=='Earlier-state span',.13,-.13)
    p<-ggplot(d,aes(estimate,y,colour=program,shape=component))+null_line()+ci_geom(aes(xmin=ci_low,xmax=ci_high))+
      geom_point(size=1.4)+facet_wrap(~dataset,nrow=1)+scale_y_continuous(breaks=1:3,labels=c('V4','V3','V2'))+
      scale_colour_manual(values=pal,guide='none')+scale_shape_manual(values=c(1,15),name=NULL)+labs(x='Aggressive-direction aligned estimate',y=NULL)
    if(v=='C') {
      z<-rd(id,'step_minus_span');z$row<-as.character(z$program)
      p2<-base_intervals(z,xlab='Step-minus-span contrast (95% CI)',facets=FALSE)+facet_wrap(~dataset,nrow=1)+theme(strip.text=element_blank(),strip.background=element_blank())
      return(p/p2+plot_layout(heights=c(1.4,1)))
    }
    p
  }
}
paired_plot <- function(id,v) {
  d<-rd(id,'paired');s<-rd(id,'summary');d$program<-factor(d$program,levels=c('V2','V3','V4'))
  if(v=='A') {
    lim<-range(c(d$first,d$second),na.rm=TRUE)
    p<-ggplot(d,aes(first,second,colour=program))+geom_abline(slope=1,intercept=0,linetype='dashed',colour=helper_colour,linewidth=.2)+
      geom_point(size=.9,alpha=.85)+facet_wrap(~program,nrow=1)+scale_colour_manual(values=pal,guide='none')+coord_cartesian(xlim=lim,ylim=lim)+
      labs(x='Earlier-state patient mean',y='LUAD patient mean')
  } else {
    long<-rbind(data.frame(program=d$program,id=d$id,state='Earlier state',value=d$first),data.frame(program=d$program,id=d$id,state='LUAD',value=d$second))
    long$state<-factor(long$state,levels=c('Earlier state','LUAD'))
    p<-ggplot(long,aes(state,value,group=id))+geom_line(colour='#A6A6A6',alpha=.45,linewidth=.17)+
      geom_point(aes(colour=program),size=.85,alpha=.85)+facet_wrap(~program,nrow=1)+scale_colour_manual(values=pal,guide='none')+labs(x=NULL,y='Patient-mean score')
  }
  d$y<-0; s$y<-0
  delta<-ggplot(d,aes(delta,y,colour=program))+geom_vline(xintercept=0,linetype='dashed',colour=helper_colour,linewidth=.2)+
    geom_point(size=.7,alpha=.8,position=position_nudge(y=.10))+ci_geom(data=s,aes(x=estimate,y=y,xmin=ci_low,xmax=ci_high),inherit.aes=FALSE,colour='#242424')+
    geom_point(data=s,aes(x=estimate,y=y),inherit.aes=FALSE,shape=18,size=1.85)+facet_wrap(~program,nrow=1)+
    scale_colour_manual(values=pal,guide='none')+scale_y_continuous(limits=c(-.15,.35),breaks=NULL)+labs(x='LUAD minus earlier states',y=NULL)+theme(axis.line.y=element_blank(),axis.ticks.y=element_blank(),strip.text=element_blank())
  p/delta+plot_layout(heights=c(3.2,1))
}
proteomics_plot <- function(v) {
  d<-rd('F4d'); d$row<-factor(d$program,levels=c('V4','V3','V2','V1'))
  effect<-ggplot(d,aes(nes,row,colour=program))+null_line()+geom_segment(aes(x=0,xend=nes,yend=row),linewidth=.3)+geom_point(shape=18,size=1.85)+
    scale_colour_manual(values=pal,guide='none')+labs(x='Observed NES',y=NULL)
  coverage<-ggplot(d,aes(coverage_fraction,row,fill=program))+geom_col(width=.4)+scale_fill_manual(values=pal,guide='none')+
    scale_x_continuous(limits=c(0,1),breaks=c(0,.5,1))+labs(x='Protein coverage',y=NULL)+theme(axis.text.y=element_blank(),axis.ticks.y=element_blank(),axis.line.y=element_blank())
  stopifnot(all(d$fdr>0))
  d$loo_label<-ifelse(is.na(d$loo_total),'NA',paste0(d$loo_significant,'/',d$loo_total))
  d$coverage_label<-paste0(d$number_measured,'/',d$number_total)
  d$fdr_label<-vapply(d$fdr,function(x)if(x>=.001)sprintf('%.3f',x)else sprintf('%.1f %%*%% 10^{%d}',x/10^floor(log10(x)),floor(log10(x))),character(1))
  metrics<-ggplot(d,aes(y=row))+geom_text(aes(x=.8,label=loo_label),size=1.95)+geom_text(aes(x=2.25,label=fdr_label),size=1.85,parse=TRUE)+geom_text(aes(x=3.6,label=coverage_label),size=1.85)+
    scale_x_continuous(limits=c(.1,4.4),breaks=c(.8,2.25,3.6),labels=c('LOO support','FDR','Measured'),position='top')+labs(x=NULL,y=NULL)+theme_void(base_family='Arial',base_size=6.2)+theme(axis.text.x=element_text(colour='#161616',size=5.8))
  if(v=='A')return(effect+metrics+plot_layout(widths=c(1.6,1)))
  fdrp<-ggplot(d,aes(-log10(fdr),row,fill=program))+geom_col(width=.4)+scale_fill_manual(values=pal,guide='none')+labs(x='-log10 FDR',y=NULL)+theme(axis.text.y=element_blank(),axis.ticks.y=element_blank(),axis.line.y=element_blank())
  support<-ggplot(d,aes(loo_significant/loo_total,row,fill=program))+geom_col(width=.4,na.rm=TRUE)+geom_text(aes(label=loo_label,x=.5),size=1.855,colour='#161616')+
    scale_fill_manual(values=pal,guide='none')+scale_x_continuous(limits=c(0,1),breaks=c(0,1))+labs(x='LOO support',y=NULL)+theme(axis.text.y=element_blank(),axis.ticks.y=element_blank(),axis.line.y=element_blank())
  effect+fdrp+coverage+support+plot_layout(widths=c(1.5,1,1,1))
}
decomposition_plot <- function(v) {
  d<-rd('F5c');d$row<-factor(d$row,levels=rev(lineage_order));d$component<-factor(d$component,levels=c('Within-cell-type','Fine-cell-type composition','Covariance','Total'))
  d$label<-paste(d$row,d$program,sep=' / '); lev<-unlist(lapply(rev(lineage_order),function(x) paste(x,c('V3','V2'),sep=' / ')));d$label<-factor(d$label,levels=lev)
  if(v=='B') {
    p<-ggplot(d,aes(estimate,label,colour=component,shape=component))+null_line()+
      ci_geom(aes(xmin=ci_low,xmax=ci_high),position=position_dodge(width=.65))+
      geom_point(position=position_dodge(width=.65),size=1.2)+scale_colour_manual(values=component_ink,name=NULL)+
      scale_shape_manual(values=c(16,15,17,18),name=NULL)+labs(x='Signed additive component / total',y=NULL)
    return(p)
  }
  d$lineage<-factor(as.character(d$row),levels=lineage_order)
  comp<-d[d$component!='Total',];tot<-d[d$component=='Total',]
  stack<-ggplot(comp,aes(estimate,label,fill=component))+geom_vline(xintercept=0,linetype='dashed',colour=helper_colour,linewidth=.2)+
    geom_col(width=.63,position=position_stack(reverse=TRUE))+scale_fill_manual(values=component_pal,name=NULL,guide=guide_legend(order=1))+
    ci_geom(data=tot,aes(x=estimate,y=label,xmin=ci_low,xmax=ci_high),inherit.aes=FALSE,colour='#222222')+
    geom_point(data=tot,aes(x=estimate,y=label,shape='Total'),inherit.aes=FALSE,size=1.3)+scale_shape_manual(values=c(Total=18),name=NULL,guide=guide_legend(order=2))+labs(x='LUAD minus earlier states',y=NULL)
  own<-ggplot(comp,aes(estimate,label,colour=component,shape=component))+null_line()+ci_geom(aes(xmin=ci_low,xmax=ci_high),position=position_dodge(width=.6))+
    geom_point(size=.9,position=position_dodge(width=.6))+scale_colour_manual(values=component_ink,guide='none')+scale_shape_manual(values=c(16,15,17),guide='none')+
    labs(x='Component estimate (95% CI)',y=NULL)+theme(axis.text.y=element_blank(),axis.ticks.y=element_blank(),axis.line.y=element_blank())
  stack<-stack+facet_grid(lineage~.,scales='free_y',space='free_y',switch='y')+scale_y_discrete(labels=function(x)sub('.* / ','',x))+theme(strip.text.y.left=element_text(angle=0,size=6.1),strip.placement='outside',panel.spacing.y=unit(1.2,'mm'))
  own<-own+facet_grid(lineage~.,scales='free_y',space='free_y')+theme(strip.text.y=element_blank(),panel.spacing.y=unit(1.2,'mm'))
  sb<-ggplot_build(stack);ob<-ggplot_build(own)
  stopifnot(nrow(sb$data[[2]])==24,nrow(sb$data[[3]])==8,nrow(ob$data[[2]])==24,nrow(ob$data[[3]])==24)
  stopifnot(max(abs(sort(sb$data[[3]]$xmin)-sort(tot$ci_low)))<1e-12,max(abs(sort(sb$data[[3]]$xmax)-sort(tot$ci_high)))<1e-12)
  stopifnot(max(abs(sort(ob$data[[2]]$xmin)-sort(comp$ci_low)))<1e-12,max(abs(sort(ob$data[[2]]$xmax)-sort(comp$ci_high)))<1e-12)
  write.csv(data.frame(panel_id='F5c',component_bars=24,total_intervals=8,own_component_intervals=24,endpoints_preserved=TRUE),file.path(outdir,'QC/F5c_GROUPED_MARK_CHECK.csv'),row.names=FALSE)
  stack+own+plot_layout(widths=c(2.5,1))
}
independent_matrix <- function(d) {
  d$row<-factor(d$row,levels=rev(lineage_order));d$col<-match(as.character(d$program),c('V2','V3','V4'));d$y<-as.numeric(d$row)
  rng<-range(c(d$ci_low,d$ci_high),na.rm=TRUE);span<-diff(rng)
  scale_local<-function(x)(x-mean(rng))/span*.77
  d$lo<-d$col+scale_local(d$ci_low);d$hi<-d$col+scale_local(d$ci_high);d$xx<-d$col+scale_local(d$estimate)
  write.csv(data.frame(range_low=rng[1],range_high=rng[2],cell_width_fraction=.77,rows=nrow(d)),file.path(outdir,'QC','F5e_INSET_SCALE.csv'),row.names=FALSE)
  ggplot(d)+geom_tile(aes(x=col,y=y,fill=estimate),width=.94,height=.82,alpha=if(publication)1 else .35)+
    geom_segment(aes(x=col+scale_local(0),xend=col+scale_local(0),y=y-.23,yend=y+.23),colour='#999999',linewidth=.22)+
    geom_segment(aes(x=lo,xend=hi,y=y,yend=y),colour='#222222',linewidth=.35)+geom_point(aes(x=xx,y=y),size=1.25)+
    scale_fill_gradient2(low=if(publication)'#D9EAF3'else'#7B5B9B',mid=if(publication)'#F9FAFE'else'white',high=if(publication)'#F0D3D3'else'#3A8790',midpoint=0,name='Difference')+
    scale_x_continuous(breaks=1:3,labels=c('V2','V3','V4'),limits=c(.45,3.55))+
    scale_y_continuous(breaks=1:4,labels=levels(d$row),limits=c(.5,4.5))+labs(x=NULL,y=NULL)+
    theme(axis.line=element_blank(),axis.ticks=element_blank(),legend.position='right',legend.key.height=unit(4,'mm'))
}

build_core <- function(s,v) {
 id<-s$panel_id;kind<-s$visualization
 if(kind=='compact_architecture') {
   d<-rd('F1e','components');z<-rd('F1e','step_minus_span');z$component<-'Step-minus-span contrast'
   d$base<-match(as.character(d$program),c('V4','V3','V2'))*.85+3.45;d$y<-d$base+ifelse(d$component=='Earlier-state span',.13,-.13)
   z$y<-match(as.character(z$program),c('V4','V3','V2'))*.85+.15
   u<-rbind(d[,c('dataset','program','component','estimate','ci_low','ci_high','y')],z[,c('dataset','program','component','estimate','ci_low','ci_high','y')])
   u$component<-factor(u$component,levels=c('Earlier-state span','Invasive-state increment','Step-minus-span contrast'))
   return(ggplot(u,aes(estimate,y,colour=program,shape=component))+null_line()+geom_hline(yintercept=3.5,colour=helper_colour,linewidth=.16)+summary_connector(id,d,'component','Earlier-state span','Invasive-state increment')+
     ci_geom(aes(xmin=ci_low,xmax=ci_high),width=.06)+geom_point(size=1.25)+facet_wrap(~dataset,nrow=1)+
     scale_colour_manual(values=pal,guide='none')+scale_shape_manual(values=c(1,15,16),name=NULL)+
     scale_y_continuous(breaks=c(6,5.15,4.3,2.7,1.85,1),labels=rep(c('V2','V3','V4'),2),limits=c(.55,6.55))+
     labs(x='Aggressive-direction estimate (95% CI)',y=NULL)+theme(legend.text=element_text(size=5.6),plot.margin=margin(.5,2,.5,1,'mm')))
 }
 if(kind=='integrated_adjusted') {
   d<-rd('F3a','integrated');d$y<-match(as.character(d$program),c('V4','V3','V2'))+ifelse(d$series=='Unadjusted',.13,-.13)
   return(ggplot(d,aes(estimate,y,colour=program,shape=series))+null_line()+summary_connector(id,d,'series','Unadjusted','V1-adjusted')+ci_geom(aes(xmin=ci_low,xmax=ci_high))+geom_point(size=1.3,stroke=.35)+
     scale_colour_manual(values=pal,guide='none')+scale_shape_manual(values=c('Unadjusted'=1,'V1-adjusted'=15),name=NULL)+facet_wrap(~dataset,nrow=1)+
     scale_y_continuous(breaks=1:3,labels=c('V4','V3','V2'))+scale_x_continuous(limits=c(-1.05,1.25),breaks=c(-1,0,1))+
     labs(x=expression('Invasive-state contrast ('*beta*')'),y=NULL)+theme(legend.position='bottom'))
 }
 if(kind=='selected_coefficients') {
   d<-rd(id);d$program<-factor(d$program,levels=c('V4','V3','V2'))
   return(ggplot(d,aes(attenuation_percent,program,colour=program))+geom_vline(xintercept=0,colour=helper_colour,linewidth=.2,linetype='dashed')+
     geom_segment(aes(x=0,xend=attenuation_percent,yend=program),linewidth=.35)+geom_point(size=1.5)+
     geom_text(aes(label=sprintf('%.1f%%',attenuation_percent)),nudge_y=.20,size=1.9,show.legend=FALSE)+
     facet_wrap(~dataset,nrow=1)+scale_colour_manual(values=pal,guide='none')+
     labs(x='Reduction in absolute coefficient after V1 adjustment (%)',y=NULL))
 }
 make_s4d<-function() {
   d<-rd('S4d');key<-c('epithelial','fibroblast_stromal','endothelial','myeloid')
   d$lineage_index<-match(d$lineage,key);stopifnot(nrow(d)==16,!anyNA(d$lineage_index),all(d$bootstrap_replicates==5000))
   centers<-c(46,35.5,25,14.5)
   d$row_y<-centers[d$lineage_index]+ifelse(d$program=='V2',2,-2)
   d$y_mm<-d$row_y+ifelse(d$component=='Within-cell-type',.58,-.58)
   domain<-range(c(d$ci_low,d$ci_high,0));padding<-diff(domain)*.055;limits<-domain+c(-padding,padding)
   breaks<-pretty(limits,n=5);breaks<-breaks[breaks>=limits[1]&breaks<=limits[2]]
   x0<-46;x1<-180;xm<-function(x)x0+(x-limits[1])/diff(limits)*(x1-x0)
   ink<-c('Within-cell-type'='#748791','Fine-cell-type composition'='#9FAAB3')
   d$x_mm<-xm(d$estimate);d$ci_low_mm<-xm(d$ci_low);d$ci_high_mm<-xm(d$ci_high)
   txt<-function(label,x,y,size=7,just='left')textGrob(label,x=unit(x,'mm'),y=unit(y,'mm'),just=just,gp=gpar(fontfamily='Arial',fontsize=size,col='#000000'))
   grobTree(
    segmentsGrob(x0=unit(xm(0),'mm'),x1=unit(xm(0),'mm'),y0=unit(10.4,'mm'),y1=unit(50.3,'mm'),gp=gpar(col='#C4CAD0',lwd=.18*72/25.4,lty='dashed')),
    txt(c('Epithelial','Fibroblast/stromal','Endothelial','Myeloid'),1,centers,7.3),
    segmentsGrob(x0=unit(30.6,'mm'),x1=unit(30.6,'mm'),y0=unit(centers-2.8,'mm'),y1=unit(centers+2.8,'mm'),gp=gpar(col='#BBC1C7',lwd=.4)),
    segmentsGrob(x0=unit(30.6,'mm'),x1=unit(31.7,'mm'),y0=unit(rep(centers,each=2)+rep(c(-2.8,2.8),4),'mm'),y1=unit(rep(centers,each=2)+rep(c(-2.8,2.8),4),'mm'),gp=gpar(col='#BBC1C7',lwd=.4)),
    txt(rep(c('V2','V3'),4),35,rep(centers,each=2)+rep(c(2,-2),4),6.7),
    segmentsGrob(x0=unit(d$ci_low_mm,'mm'),x1=unit(d$ci_high_mm,'mm'),y0=unit(d$y_mm,'mm'),y1=unit(d$y_mm,'mm'),gp=gpar(col=ink[d$component],lwd=.36*72/25.4,lineend='butt')),
    segmentsGrob(x0=unit(rep(c(d$ci_low_mm,d$ci_high_mm),each=1),'mm'),x1=unit(c(d$ci_low_mm,d$ci_high_mm),'mm'),y0=unit(rep(d$y_mm,2)-.38,'mm'),y1=unit(rep(d$y_mm,2)+.38,'mm'),gp=gpar(col=rep(ink[d$component],2),lwd=.36*72/25.4)),
    pointsGrob(x=unit(d$x_mm,'mm'),y=unit(d$y_mm,'mm'),pch=ifelse(d$component=='Within-cell-type',16,15),size=unit(1.5,'mm'),gp=gpar(col=ink[d$component],fill=ink[d$component])),
    segmentsGrob(x0=unit(x0,'mm'),x1=unit(x1,'mm'),y0=unit(8,'mm'),y1=unit(8,'mm'),gp=gpar(col='#000000',lwd=.30*72/25.4)),
    segmentsGrob(x0=unit(xm(breaks),'mm'),x1=unit(xm(breaks),'mm'),y0=unit(8,'mm'),y1=unit(7,'mm'),gp=gpar(col='#000000',lwd=.26*72/25.4)),
    txt(format(breaks,trim=TRUE),xm(breaks),5.8,6.7,'centre'),
    txt('V1-adjusted component estimate (95% CI)',mean(c(x0,x1)),1.9,7.4,'centre'),
    pointsGrob(x=unit(c(61,106),'mm'),y=unit(rep(55.8,2),'mm'),pch=c(16,15),size=unit(1.5,'mm'),gp=gpar(col=ink,fill=ink)),
    txt('Within-cell-type',64,55.8,6.5),txt('Fine-cell-type composition',109,55.8,6.5))
 }
 if(kind=='v1_adjusted_decomposition')return(make_s4d())
 if(kind=='compact_vi')return(base_intervals(rd('F3b','vi'),'Microscopic vascular invasion\nOR (95% CI)',log=TRUE))
 if(kind=='proteomics_strip') {
   d<-rd('F4d');d$row<-factor(d$program,levels=c('V4','V3','V2','V1'));d$loo<-ifelse(is.na(d$loo_total),'NA',paste0(d$loo_significant,'/',d$loo_total))
   p<-ggplot(d,aes(nes,row,colour=program))+null_line()+geom_segment(aes(x=0,xend=nes,yend=row),linewidth=.3)+geom_point(shape=18,size=1.6)+
     scale_colour_manual(values=pal,guide='none')+labs(x='Observed NES',y=NULL)
   q<-ggplot(d,aes(.5,row))+geom_text(aes(label=loo),size=1.95)+scale_x_continuous(limits=c(0,1),breaks=.5,labels='LOO support',position='top')+
     labs(x=NULL,y=NULL)+theme_void(base_family='Arial',base_size=6.2)+theme(axis.text.x=element_text(size=5.6,colour='#161616'))
   return(p+q+plot_layout(widths=c(3,1)))
 }
 if(kind=='pooled_robustness') {
   d<-rd('S3c','robustness');d$row<-factor(d$row,levels=rev(c('Primary (3 cohorts)','Modified HK sensitivity','Omit GSE50081','Omit GSE31210','Omit GSE30219')))
   stopifnot(all(d$ci_low>0))
   return(ggplot(d,aes(estimate,row,colour=program,shape=series))+null_line(TRUE)+ci_geom(aes(xmin=ci_low,xmax=ci_high))+geom_point(size=1.2)+
     facet_wrap(~program,nrow=1)+scale_colour_manual(values=pal,guide='none')+
     scale_shape_manual(values=c('Ordinary HK'=18,'Modified HK'=1,'Leave-one-cohort-out'=16),labels=c('Ordinary HK'='Primary','Modified HK'='Modified HK sensitivity','Leave-one-cohort-out'='Cohort omission'),name=NULL)+
     scale_x_log10(breaks=c(.1,1,10),labels=c('0.1','1','10'))+labs(x='OS hazard ratio (95% CI)',y=NULL)+theme(legend.position='bottom',axis.text.y=element_text(size=5.5)))
 }
 if(kind=='coverage') {
   d<-rd(id);d$row<-factor(d$row,levels=rev(unique(d$row)));d$label<-paste0(sprintf('%.0f%%',d$coverage*100),'\n',d$mapped,'/',d$total)
   return(ggplot(d,aes(program,row))+geom_tile(fill=if(publication)'#F9FAFE'else'#F3F4F6',colour=if(publication)'#E7E9EF'else'white',linewidth=.4)+geom_text(aes(label=label),size=1.85)+
     labs(x=NULL,y=NULL)+theme(axis.line=element_blank(),axis.ticks=element_blank()))
 }

 if(kind=='schematic') {
   svg<-file.path(root,'data_processed','figure_assets','F1a.svg')
   stopifnot(file.exists(svg))
   pngfile<-tempfile(fileext='.png')
   rsvg::rsvg_png(svg,file=pngfile,width=1400,height=1400)
   im<-png::readPNG(pngfile);unlink(pngfile)
   return(rasterGrob(im,width=unit(62,'mm'),height=unit(62,'mm'),interpolate=TRUE))
 }
 if(kind=='trajectory') {
   d<-rd(id,'raw');z<-rd(id,'summary');levels<-c('Normal','AAH','AIS','MIA','LUAD');d$x<-match(d$stage,levels);z$x<-match(z$stage,levels)
   # Geometric x-only separation leaves every measured score unchanged.
   d<-d[order(d$program,d$stage,d$sample),];ids<-sort(unique(d$patient));d$offset<-(match(d$patient,ids)-(length(ids)+1)/2)/length(ids)*.24
   return(ggplot(d,aes(x+offset,value,colour=program))+geom_hline(yintercept=0,colour=helper_colour,linewidth=.2)+geom_point(size=.5,alpha=.35)+
     geom_line(data=z,aes(x=x,y=estimate,group=program),linewidth=.2,alpha=.5)+geom_errorbar(data=z,aes(x=x,y=estimate,ymin=ci_low,ymax=ci_high),width=.12,linewidth=.3)+
     geom_point(data=z,aes(x=x,y=estimate),shape=15,size=1.3)+facet_wrap(~program,nrow=1)+scale_colour_manual(values=pal,guide='none')+
     scale_x_continuous(breaks=1:5,labels=levels)+labs(x='Pathology category',y='Programme score')+theme(axis.text.x=element_text(size=5.7)))
 }
 if(kind=='models')return(models_plot(id,v))
 if(kind=='null') {
   d<-rd(id,'null');z<-rd(id,'observed')
   stopifnot(all(z$p_emp>0));z$label<-vapply(z$p_emp,function(x)sprintf('italic(P) == %.1f %%*%% 10^{%d}',x/10^floor(log10(x)),floor(log10(x))),character(1))
   return(ggplot(d,aes(null_stat,fill=program))+geom_histogram(bins=55,linewidth=0,alpha=.24)+geom_vline(data=z,aes(xintercept=observed_stat,colour=program),linewidth=.65,linetype='dashed')+
     geom_text(data=z,aes(x=-Inf,y=Inf,label=label),inherit.aes=FALSE,parse=TRUE,hjust=-.03,vjust=1.3,size=1.85)+facet_wrap(~program,nrow=1,scales='free_y')+
     scale_fill_manual(values=pal,guide='none')+scale_colour_manual(values=pal,guide='none')+labs(x='AICc gradual minus invasive step',y='Matched sets'))
 }
 if(kind=='architecture_components')return(component_plot(id,v))
 if(kind=='populations') {
   d<-rd(id,'counts')
   if(v=='A') {
     d$population<-factor(d$population,levels=rev(d$population))
     return(ggplot(d,aes(n,population))+geom_col(fill=if(publication)'#8DC4D9'else'#617F8C',width=.5)+geom_text(aes(label=n),hjust=-.2,size=2.2)+scale_x_continuous(limits=c(0,45),breaks=c(0,20,40))+labs(x='Patients',y=NULL))
   }
   m<-rd(id,'membership');cols<-c('Direct contrast','Stage IA','Paired');names(m)<-gsub('\\.',' ',names(m))
   key<-apply(m[,cols,drop=FALSE],1,paste0,collapse='');tab<-sort(table(key),decreasing=TRUE);keys<-names(tab)
   d2<-data.frame(key=factor(keys,levels=keys),n=as.numeric(tab))
   p1<-ggplot(d2,aes(key,n))+geom_col(fill=if(publication)'#8DC4D9'else'#617F8C',width=.65)+geom_text(aes(label=n),vjust=-.2,size=2)+scale_y_continuous(expand=expansion(mult=c(0,.18)))+labs(x=NULL,y='Patients')+theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())
   m2<-expand.grid(key=keys,subset=cols,stringsAsFactors=FALSE);m2$value<-mapply(function(k,j)as.integer(substr(k,j,j)),m2$key,match(m2$subset,cols));m2$key<-factor(m2$key,levels=keys);m2$subset<-factor(m2$subset,levels=rev(cols))
   p2<-ggplot(m2,aes(key,subset))+geom_point(colour='#DDDDDD',size=1.5)+geom_point(data=m2[m2$value==1,],colour='#2A5366',size=1.5)+labs(x=NULL,y=NULL)+theme_void(base_family='Arial',base_size=6)+theme(axis.text.y=element_text(size=5.7))
   return((p1+labs(subtitle='40 patients overall')+theme(plot.subtitle=element_text(size=6))) / p2+plot_layout(heights=c(2,1)))
 }
 if(kind=='paired')return(paired_plot(id,v))
 if(kind=='adjusted') {
   d<-rd(id);d$y<-match(as.character(d$program),c('V4','V3','V2'))+ifelse(d$series=='Unadjusted',.12,-.12)
   return(ggplot(d,aes(estimate,y,colour=program,shape=series))+null_line()+ci_geom(aes(xmin=ci_low,xmax=ci_high))+geom_point(size=1.3,stroke=.35)+
     scale_colour_manual(values=pal,guide='none')+scale_shape_manual(values=c('Unadjusted'=1,'V1-adjusted'=15),name=NULL)+
     scale_y_continuous(breaks=1:3,labels=c('V4','V3','V2'))+scale_x_continuous(limits=c(-1.05,1.25),breaks=c(-1,0,1))+
     labs(x=expression('Invasive-state contrast ('*beta*')'),y=NULL)+theme(legend.position='none'))
 }
 if(kind=='log_interval')return(base_intervals(rd(id),'Odds ratio per 1-SD increase',log=TRUE))
 if(kind=='interval')return(base_intervals(rd(id),switch(id,F2b=expression('Invasive-state contrast ('*beta*')'),F2e=expression('Stage IA invasive-state contrast ('*beta*')'),F4b=expression('Solid minus pure non-solid ('*beta*')'),F5d=expression('Invasive-state contrast ('*beta*')'),S2b='LUAD minus AAH',S5c='IAC minus pooled AIS/MIA','Contrast')))
 if(kind=='interval_matrix') {
   d<-rd(id);lg<-id=='S3b'
   if(id=='S3b') {
     stopifnot(all(d$ci_low>0));d$row<-factor(d$row,levels=rev(unique(d$row)));d$y<-as.numeric(d$row)+(match(d$program,c('V2','V3','V4'))-2)*.18
     return(ggplot(d,aes(estimate,y,colour=program,shape=program))+null_line(TRUE)+ci_geom(aes(xmin=ci_low,xmax=ci_high))+geom_point(size=1.25)+
       scale_colour_manual(values=pal,name=NULL)+scale_shape_manual(values=c(V2=16,V3=15,V4=17),name=NULL)+
       scale_x_log10(breaks=c(.1,1,10,100),labels=c('0.1','1','10','100'))+
       scale_y_continuous(breaks=seq_along(levels(d$row)),labels=gsub(': ',':\n',levels(d$row)))+
       labs(x='OS hazard ratio (95% CI)',y=NULL)+theme(axis.text.y=element_text(size=5.4)))
   }
   if(lg)d$row<-gsub('V1 ADJUSTED PURE','V1 adjusted',gsub('V1 ADJUSTED','V1 adjusted',gsub('CLINICALLY ADJUSTED','Clinically adjusted',d$row)))
   axis_title<-switch(id,F4a='Invasive minus AIS/MIA (standardized score)',F5b=expression('V1-adjusted invasive-state contrast ('*beta*')'),S2c='Cohort-specific standardized contrast',S3b='Overall-survival hazard ratio (95% CI)','Effect (95% CI)')
   return(base_intervals(d,axis_title,facets=TRUE,log=lg))
 }
 if(kind=='regional_pairs') {
   d<-rd(id);d$row<-factor(d$id,levels=rev(unique(d$id)))
   return(ggplot(d,aes(y=row,colour=program))+geom_segment(aes(x=first,xend=second,yend=row),linewidth=.3,alpha=.55)+
     geom_point(aes(x=first),shape=1,size=1.35,stroke=.35)+geom_point(aes(x=second),shape=16,size=1.35)+facet_wrap(~program,nrow=1)+
     scale_colour_manual(values=pal,guide='none')+labs(x='Regional program score',y=NULL))
 }
 if(kind=='proteomics')return(proteomics_plot(v))
 if(kind=='survival') {
   d<-rd(id);d<-order_row(d)
   if('pooled'%in%names(d)) {
     d$shape<-ifelse(d$pooled,'Pooled','Cohort')
     p<-ggplot(d,aes(estimate,row,colour=program))+null_line(TRUE)+ci_geom(aes(xmin=ci_low,xmax=ci_high))+geom_point(aes(shape=shape),size=1.5)+
       scale_shape_manual(values=c(Cohort=16,Pooled=18),guide='none')+scale_colour_manual(values=pal,guide='none')+facet_wrap(~program,nrow=1)+
       scale_x_log10(breaks=c(.03,.1,.3,1,3,10,30))+labs(x=if(id=='S3d')'Recurrence-free survival HR (95% CI)' else 'Overall-survival hazard ratio (95% CI)',y=NULL)
     if(FALSE && id=='F4e') {
       # n/events labels are derived from unique cohort rows, with a deterministic pooled sum.
       labelmap<-setNames(paste0(as.character(d$row),'  (',d$n,'/',d$events,')'),as.character(d$row))
       p<-p+scale_y_discrete(labels=labelmap)
     }
     return(p)
   }
   return(base_intervals(d,'OS hazard ratio after cohort omission',facets=TRUE,log=TRUE))
 }
 if(kind=='umap') {
   cache<-file.path(datadir,'F5a_umap_point_layer.png')
   stopifnot(file.exists(cache))
   if(!file.exists(cache)) {
     d<-rd('F5a','umap'); cols<-c('Epithelial'='#D58E69','Fibroblast/stromal'='#C4AB54','Endothelial'='#6294A3','Myeloid'='#9B7FA7','T / NK'='#80AD9C','B / plasma'='#B1B86F','Mast cells'='#B56F70','Other'='#B7B7B7')
     grDevices::png(cache,width=1900,height=1350,res=600,type='cairo',bg='white')
     par(mar=c(0,0,0,0),xaxs='i',yaxs='i');plot(range(d$UMAP1),range(d$UMAP2),type='n',axes=FALSE,xlab='',ylab='',asp=1)
     for(g in names(cols)){ii<-d$display_lineage==g;points(d$UMAP1[ii],d$UMAP2[ii],pch=16,cex=.085,col=adjustcolor(cols[g],alpha.f=.43))}
     grDevices::dev.off()
     write.csv(data.frame(nuclei=nrow(d),display_categories=length(unique(d$display_lineage)),coordinate_transform='NONE',point_exclusion=0),file.path(root,'results/cellular','UMAP_RENDER_COUNTS.csv'),row.names=FALSE)
   }
   im<-png::readPNG(cache)
   labs<-data.frame(x=c(.66,.78,.44,.34,.49,.73,.10),y=c(.87,.55,.70,.40,.10,.26,.66),label=c('Epithelial','Fibroblast/\nstromal','Endothelial','Myeloid','T / NK','B / plasma','Mast cells'))
   image_h<-41;image_w<-image_h*dim(im)[2]/dim(im)[1];frac<-image_w/82
   labs$x<-.5+(labs$x-.5)*frac
   g<-grobTree(rasterGrob(im,width=unit(image_w,'mm'),height=unit(image_h,'mm'),interpolate=TRUE),
     segmentsGrob(x0=.05,x1=.15,y0=.045,y1=.045,gp=gpar(col='#555555',lwd=.4)),segmentsGrob(x0=.05,x1=.05,y0=.045,y1=.16,gp=gpar(col='#555555',lwd=.4)),
     textGrob('UMAP1',x=.16,y=.045,just='left',gp=gpar(fontfamily='Arial',fontsize=5.2)),textGrob('UMAP2',x=.05,y=.18,just='bottom',gp=gpar(fontfamily='Arial',fontsize=5.2)),
     textGrob(labs$label,x=labs$x,y=labs$y,gp=gpar(fontfamily='Arial',fontsize=5.6,col='#161616')))
   return(g)
 }
 if(kind=='decomposition')return(decomposition_plot(v))
 if(kind=='independent') {
   d<-rd(id)
   if(v=='B')return(independent_matrix(d))
   return(base_intervals(d,'IAC minus pooled AIS/MIA (95% CI)',facets=TRUE))
 }
 if(kind=='heldout') {
   d<-rd(id);d$model<-factor(d$model,levels=c('Gradual','Invasive step','Hinge'))
   z<-aggregate(value~program+model,d,mean)
   return(ggplot(d,aes(model,value,group=patient))+geom_line(colour='#BDBDBD',alpha=.23,linewidth=.13)+geom_point(aes(colour=program),size=.55,alpha=.7)+
     geom_point(data=z,aes(model,value,colour=program),inherit.aes=FALSE,shape=15,size=1.5)+facet_wrap(~program,nrow=1)+scale_colour_manual(values=pal,guide='none')+labs(x='Model',y='Predictive log-likelihood')+theme(axis.text.x=element_text(size=5.6)))
 }
 if(kind=='scenario_intervals') {
   d<-rd(id);d$row<-as.character(d$program);d<-order_row(d,c('V2','V3','V4'));d$y<-as.numeric(d$row)+ifelse(d$series=='Full sequence',.12,-.12)
   return(ggplot(d,aes(estimate,y,colour=program,shape=series))+null_line()+ci_geom(aes(xmin=ci_low,xmax=ci_high))+geom_point(size=1.3)+facet_wrap(~dataset,nrow=1)+scale_colour_manual(values=pal,guide='none')+
     scale_shape_manual(values=c('Full sequence'=16,'Lesion-only'=1),name=NULL)+scale_y_continuous(breaks=1:3,labels=c('V4','V3','V2'))+labs(x='Aggressive-direction step-minus-span contrast (95% CI)',y=NULL))
 }
 if(kind=='patient_influence') {
   d<-rd(id);d$scenario<-factor(d$scenario,levels=c('Full sequence','Lesion-only'))
   d$offset<-ave(d$estimate,interaction(d$dataset,d$scenario,d$program),FUN=function(x)seq(-.08,.08,length.out=length(x)))
   d$y<-match(as.character(d$program),c('V4','V3','V2'))+d$offset
   d$block<-factor(paste(d$scenario,d$dataset,sep=' | '),levels=c('Full sequence | GSE308103','Full sequence | GSE319666','Lesion-only | GSE308103','Lesion-only | GSE319666'))
   return(ggplot(d,aes(estimate,y,colour=program))+null_line()+geom_point(size=.6,alpha=.6)+
     facet_wrap(~block,ncol=2)+scale_colour_manual(values=pal,guide='none')+theme(strip.text=element_text(size=5.8,angle=0))+scale_y_continuous(breaks=1:3,labels=c('V4','V3','V2'))+labs(x='Leave-one-patient-out step-minus-span contrast',y=NULL))
 }
 if(kind=='alternative') {
   d<-rd(id);d$row<-as.character(d$program)
   return(base_intervals(d,'Alternative-score step-minus-span contrast (95% CI)')+geom_point(aes(x=primary,shape='Primary reference'),size=1.6)+scale_shape_manual(values=c('Primary reference'=1),name=NULL)+facet_wrap(~scenario,nrow=1)+theme(legend.position='none'))
 }
 if(kind=='stage_cross') {
   d<-rd(id);d$row<-gsub('minus','minus',d$row)
   d$facet<-gsub('*','',d$facet,fixed=TRUE)
   return(base_intervals(d,'Difference in standardized score',facets=FALSE)+geom_text(data=d[!is.finite(d$estimate),],aes(x=0,label='NE'),size=1.9,colour='#555555')+facet_grid(facet~program)+theme(axis.text.y=element_text(size=5.5)))
 }
 if(kind=='coverage') {
   d<-rd(id);d$row<-factor(d$row,levels=rev(unique(d$row)));d$label<-sprintf('%.0f%%',d$coverage*100)
   if(id=='S3a')d$label<-paste0(d$label,'\n',d$mapped,'/',d$total)
   p<-ggplot(d,aes(program,row,fill=coverage))+geom_tile(colour='white',linewidth=.4)+geom_text(aes(label=label),size=if(id=='S3a')1.85 else 2)+
     scale_fill_gradient(low='#F2F1EC',high='#4E7585',limits=c(0,1),name='Coverage')+labs(x=NULL,y=NULL)+
     theme(axis.line=element_blank(),axis.ticks=element_blank(),legend.position='bottom',legend.key.width=unit(10,'mm'))
   if(id=='S3a')p<-p+geom_tile(data=d[d$coverage<.8,],fill=NA,colour='#202020',linewidth=.45)+labs(caption='Outlined cells: below 80% member coverage')
   return(p)
 }
 if(kind=='cindex') {
   d<-rd(id);d$model<-factor(d$model,levels=c('A','B','C','D'),labels=c('Clinical','Clinical +\nORACLE','Clinical +\nORACLE + V2','Clinical +\nORACLE + V3'))
   return(ggplot(d,aes(model,value,group=row,colour=row))+geom_line(linewidth=.3)+geom_point(size=1.3)+scale_colour_brewer(palette='Dark2',name=NULL)+labs(x=NULL,y='Apparent C-index')+theme(axis.text.x=element_text(size=5.4))+theme(legend.position='bottom'))
 }
 if(kind=='hk') {
   d<-rd(id);d$row<-as.character(d$program);d<-order_row(d,c('V2','V3','V4'));d$y<-as.numeric(d$row)+ifelse(d$series=='Ordinary HK',.12,-.12)
   return(ggplot(d,aes(estimate,y,colour=program,shape=series))+null_line(TRUE)+ci_geom(aes(xmin=ci_low,xmax=ci_high))+geom_point(size=1.3)+scale_x_log10()+scale_colour_manual(values=pal,guide='none')+
     scale_shape_manual(values=c('Ordinary HK'=16,'Modified HK'=1),labels=c('Ordinary HK'='Primary HK','Modified HK'='Modified HK (sensitivity)'),name=NULL)+scale_y_continuous(breaks=1:3,labels=c('V4','V3','V2'))+labs(x='Pooled OS HR (95% CI)',y=NULL)+theme(legend.position='bottom'))
 }
 if(kind=='eligibility') {
   d<-rd(id)
   if(v=='B')return(ggplot(d,aes(pre,luad,shape=eligible,colour=lineage))+geom_hline(yintercept=6,linetype='dashed',colour=helper_colour,linewidth=.25)+geom_vline(xintercept=6,linetype='dashed',colour=helper_colour,linewidth=.25)+
       geom_point(size=1.4,alpha=.7)+facet_wrap(~lineage,nrow=1)+scale_colour_brewer(palette='Set2',guide='none')+scale_shape_manual(values=c('FALSE'=1,'TRUE'=16),labels=c('FALSE'='No','TRUE'='Yes'),name='Eligible')+labs(x='Earlier-state patients',y='LUAD patients'))
   long<-rbind(data.frame(row=d$row,lineage=d$lineage,group='Earlier state',n=d$pre),data.frame(row=d$row,lineage=d$lineage,group='LUAD',n=d$luad))
   long$group<-factor(long$group,levels=c('Earlier state','LUAD'))
   # Four lineage blocks preserve all 43 fine-type labels at native size.
   long$row<-factor(long$row,levels=rev(unique(d$row)))
   return(ggplot(long,aes(group,row,fill=n))+geom_tile(colour='white',linewidth=.2)+geom_text(aes(label=n),size=1.85)+facet_wrap(~lineage,nrow=1,scales='free_y')+
     scale_fill_gradient(low='#F3F3EF',high='#567B8B',name='Patients')+labs(x=NULL,y=NULL,caption='Fine-type eligibility: at least 6 patients in each group')+theme(axis.line=element_blank(),axis.ticks=element_blank(),axis.text.y=element_text(size=5.5),legend.position='bottom'))
 }
 if(kind=='sensitivity_grid') {
   d<-rd(id);d$condition[d$condition=='Without TMM normalization']<-'Without TMM\nnormalization';d$condition[d$condition=='Alternative fine-type scores']<-'Alternative\nfine-type scores';d$condition[d$condition=='Leave-one-patient-out']<-'Leave-one-\npatient-out';d$column<-paste(d$program,ifelse(d$component=='Within-cell-type','Within-cell-type','Total'));d$condition<-factor(d$condition,levels=rev(unique(d$condition)))
   return(ggplot(d,aes(column,condition,fill=value))+geom_tile(colour='white',linewidth=.5)+geom_text(aes(label=sprintf('%.3f',value)),size=1.9)+scale_fill_gradient(low=if(publication)'#F9FAFE'else'#F2F0EB',high=if(publication)'#8DC4D9'else'#507383',name='Max. abs. deviation')+
     labs(x=NULL,y=NULL)+theme(axis.line=element_blank(),axis.ticks=element_blank(),axis.text.x=element_text(angle=30,hjust=1),legend.position='bottom'))
 }
 if(kind=='decomp_lopo') {
   d<-rd(id);d<-order_row(d);d$y<-as.numeric(d$row)+ifelse(d$component=='Total',.13,-.13)
   return(ggplot(d,aes(estimate,y,colour=component,shape=component))+null_line()+ci_geom(aes(xmin=ci_low,xmax=ci_high))+geom_point(size=1.1)+facet_wrap(~program,nrow=2)+
     scale_colour_manual(values=component_ink,name=NULL)+scale_shape_manual(values=c('Total'=18,'Within-cell-type'=16),name=NULL)+scale_y_continuous(breaks=1:4,labels=rev(lineage_order))+labs(x='Leave-one-patient-out range',y=NULL))
 }
 if(kind=='independent_adjusted') {
   d<-rd(id);d<-order_row(d);d$y<-as.numeric(d$row)+ifelse(d$series=='Primary',.12,-.12)
   return(ggplot(d,aes(estimate,y,colour=program,shape=series))+null_line()+ci_geom(aes(xmin=ci_low,xmax=ci_high))+geom_point(aes(size=series))+scale_size_manual(values=c('Primary'=1.2,'V1 residualized'=.95),guide='none')+facet_wrap(~program,nrow=1)+
     scale_colour_manual(values=pal,guide='none')+scale_shape_manual(values=c('Primary'=16,'V1 residualized'=1),labels=c('Primary'='Primary','V1 residualized'='V1-residualized sensitivity'),name=NULL)+scale_y_continuous(breaks=1:4,labels=rev(lineage_order))+labs(x='Independent patient-level contrast',y=NULL))
 }
 if(kind=='concordance') {
   d<-rd(id);d$patient<-factor(d$patient,levels=rev(unique(d$patient)));d$column<-factor(d$column,levels=unique(d$column));d$label<-ifelse(d$match==1,'+','-')
   d$lineage<-factor(sub('^V[23] ','',as.character(d$column)),levels=lineage_order);d$programme<-factor(sub(' .*','',as.character(d$column)),levels=c('V2','V3'))
   p<-ggplot(d,aes(programme,patient,fill=factor(match)))+geom_tile(colour='white',linewidth=.4)+geom_text(aes(label=label),size=2.1)+
     scale_fill_manual(values=if(publication)c('0'='#F1DADA','1'='#DCECF3')else c('0'='#DFD7E6','1'='#C1D8D6'),guide='none')+labs(x=NULL,y=NULL)+theme(axis.line=element_blank(),axis.ticks=element_blank(),axis.text.x=element_text(angle=0,hjust=.5,size=5.6),strip.text=element_text(size=5.6))+facet_grid(~lineage)
   if(v=='B') {
     z<-aggregate(match~patient,d,sum)
     pb<-ggplot_build(p);stopifnot(nrow(pb$data[[1]])==72,nrow(pb$data[[2]])==72,nrow(z)==9)
     stopifnot(all(paste(d$programme,d$lineage)==as.character(d$column)))
     write.csv(data.frame(panel_id='S5b',matrix_cells=72,direction_symbols=72,patient_counts=9,identities_preserved=TRUE),file.path(outdir,'QC/S5b_GROUPED_MARK_CHECK.csv'),row.names=FALSE)
     q<-ggplot(z,aes(match,patient))+geom_col(fill=if(publication)'#8DC4D9'else'#567D89',width=.55)+geom_text(aes(label=paste0(match,'/8')),hjust=-.1,size=1.85)+scale_x_continuous(limits=c(0,9.5),breaks=c(0,4,8))+labs(x='Concordant',y=NULL)+theme(axis.text.y=element_blank(),axis.ticks.y=element_blank(),axis.line.y=element_blank())
     return(p+q+plot_layout(widths=c(4,1)))
   }
   return(p)
 }
 stop('Unknown panel kind: ',kind)
}


as_grob<-function(p){if(inherits(p,'grob'))p else if(inherits(p,'patchwork'))patchworkGrob(p)else ggplotGrob(p)}
cores<-list()
draw_panel<-function(id,x=0,y=0,w,h,page_h=h){
  ss<-spec[spec$panel_id==id,];g<-cores[[id]]
  pushViewport(viewport(x=unit(x,'mm'),y=unit(page_h-y-h,'mm'),width=unit(w,'mm'),height=unit(h,'mm'),just=c('left','bottom')))
  grid.text(ss$letter,x=unit(.5,'mm'),y=unit(h-.5,'mm'),just=c('left','top'),gp=gpar(fontfamily='Arial',fontsize=8,fontface='bold',col='#161616'))
  top<-3
  pushViewport(viewport(x=unit(0,'mm'),y=unit(0,'mm'),width=unit(w,'mm'),height=unit(h-top,'mm'),just=c('left','bottom')))
  grid.draw(g);popViewport(2)
}
save_draw<-function(draw,stem,w,h,dpi=600,svg=TRUE){
  grDevices::cairo_pdf(paste0(stem,'.pdf'),width=w/25.4,height=h/25.4,family='Arial',bg='white');draw();grDevices::dev.off()
  if(svg){svglite::svglite(paste0(stem,'.svg'),width=w/25.4,height=h/25.4,system_fonts=list(sans='Arial'),bg='white');draw();grDevices::dev.off()}
  grDevices::png(paste0(stem,'.png'),width=w/25.4,height=h/25.4,units='in',res=dpi,type='cairo',bg='white');draw();grDevices::dev.off()
}
layouts<-list(
 F1=list(size=c(183,172),panels=rbind(c('F1a',0,0,75,62),c('F1b',75,0,108,62),c('F1c',0,66,85,46),c('F1d',85,66,98,46),c('F1e',7.5,116,168,52))),
 F2=list(size=c(183,127),panels=rbind(c('F2a',0,0,49,52),c('F2b',49,0,63,52),c('F2c',112,0,71,52),c('F2d',0,56,115,71),c('F2e',115,56,68,71))),
 F3=list(size=c(183,138),panels=rbind(c('F3a',0,0,122,62),c('F3b',122,0,61,62),c('F3c',22.5,66,138,68))),
 F4=list(size=c(183,167),panels=rbind(c('F4a',0,0,100,40),c('F4b',100,0,83,40),c('F4c',0,42,110,68),c('F4d',110,42,73,68),c('F4e',0,112,183,55))),
 F5=list(size=c(183,170),panels=rbind(c('F5a',0,0,82,46),c('F5b',82,0,101,46),c('F5c',0,48,183,70),c('F5d',0,120,68,49),c('F5e',68,120,115,49))),
 S1=list(size=c(183,205),panels=rbind(c('S1a',0,0,183,42),c('S1b',7.5,46,168,48),c('S1c',7.5,98,168,49),c('S1d',10.5,151,162,50))),
 S2=list(size=c(183,116),panels=rbind(c('S2a',0,0,183,58),c('S2b',0,62,75,50),c('S2c',75,62,108,50))),
 S3=list(size=c(183,161),panels=rbind(c('S3a',0,0,88,48),c('S3b',88,0,95,48),c('S3c',0,52,183,58),c('S3d',0,114,183,43))),
 S4=list(size=c(183,214),panels=rbind(c('S4a',0,0,183,64),c('S4b',0,68,91.5,78),c('S4c',91.5,68,91.5,78),c('S4d',0,150,183,64))),
 S5=list(size=c(183,172),panels=rbind(c('S5a',0,0,183,53),c('S5b',0,57,183,53),c('S5c',20.5,114,142,54))))
coords<-do.call(rbind,lapply(names(layouts),function(f){d<-as.data.frame(layouts[[f]]$panels);names(d)<-c('panel_id','x_mm','y_mm','width_mm','height_mm');d$figure_id<-f;d}))
write.csv(coords,file.path(outdir,'LAYOUT_COORDINATES.csv'),row.names=FALSE)
write.csv(data.frame(figure_id=names(layouts),width_mm=183,height_mm=sapply(layouts,function(x)x$size[2])),file.path(outdir,'FIGURE_SIZES.csv'),row.names=FALSE)
figures<-if(length(only_fig))only_fig else names(layouts)
short_titles<-c(S1a='Held-out prediction',S1b='Lesion-only precision',S1c='Lesion-only components',S1d='Patient influence',S2a='Stage-specific effects',S2b='AAH to LUAD',S2c='Cohort-specific effects',S3a='Programme coverage',S3b='Adjusted OS',S3c='Pooled OS robustness',S3d='Recurrence-free survival',S4a='Fine-type support',S4b='Decomposition sensitivity',S4c='Patient influence',S4d='V1-adjusted decomposition',S5a='V1 sensitivity',S5b='Directional robustness',S5c='Spatial support')
spec$title[spec$panel_id%in%names(short_titles)]<-short_titles[spec$panel_id[spec$panel_id%in%names(short_titles)]]
for(fig in figures){
  cat('RENDER',phase,fig,'\n');flush.console();ss<-spec[spec$figure_id==fig,]
  for(i in seq_len(nrow(ss))){si<-ss[i,];cores[[si$panel_id]]<-as_grob(build_core(si,si$layout));l<-coords[coords$panel_id==si$panel_id,];local({id<-si$panel_id;w<-as.numeric(l$width_mm);h<-as.numeric(l$height_mm);save_draw(function(){grid.newpage();draw_panel(id,w=w,h=h)},file.path(outdir,'PANELS',id),w,h)})}
  if(fig=='F2'){
    stopifnot(length(cores[['F2b']]$heights)==length(cores[['F2c']]$heights))
    hh<-unit.pmax(cores[['F2b']]$heights,cores[['F2c']]$heights)
    cores[['F2b']]$heights<-hh;cores[['F2c']]$heights<-hh
    layout_check<-do.call(rbind,lapply(c('F2b','F2c'),function(id){g<-cores[[id]];l<-coords[coords$panel_id==id,];data.frame(panel_id=id,box_width_mm=as.numeric(l$width_mm),box_height_mm=as.numeric(l$height_mm),effective_panel_height_mm=as.numeric(l$height_mm)-3-sum(convertHeight(g$heights,'mm',valueOnly=TRUE)),effective_panel_width_mm=as.numeric(l$width_mm)-sum(convertWidth(g$widths,'mm',valueOnly=TRUE)))}))
    write.csv(layout_check,file.path(outdir,'QC/F2_PLOT_REGION_CHECK.csv'),row.names=FALSE)
    for(id in c('F2b','F2c'))local({pid<-id;l<-coords[coords$panel_id==pid,];w<-as.numeric(l$width_mm);h<-as.numeric(l$height_mm);save_draw(function(){grid.newpage();draw_panel(pid,w=w,h=h)},file.path(outdir,'PANELS',pid),w,h)})
  }
  lay<-layouts[[fig]];W<-lay$size[1];H<-lay$size[2]
  draw_whole<-function(newpage=TRUE){if(newpage)grid.newpage();for(j in seq_len(nrow(lay$panels))){l<-lay$panels[j,];draw_panel(l[1],x=as.numeric(l[2]),y=as.numeric(l[3]),w=as.numeric(l[4]),h=as.numeric(l[5]),page_h=H)}}
  save_draw(function()draw_whole(),file.path(outdir,'FIGURES',fig),W,H)
  grDevices::cairo_pdf(file.path(outdir,'ACTUAL_SIZE_PROOFS',paste0(fig,'_ACTUAL_SIZE.pdf')),width=210/25.4,height=297/25.4,family='Arial')
  grid.newpage();grid.text(paste(fig,'| print at 100%; no fit-to-page'),x=unit(13.5,'mm'),y=unit(290,'mm'),just=c('left','top'),gp=gpar(fontfamily='Arial',fontsize=8))
  pushViewport(viewport(x=unit(13.5,'mm'),y=unit(297-18-H,'mm'),width=unit(W,'mm'),height=unit(H,'mm'),just=c('left','bottom')));draw_whole(FALSE);popViewport()
  grid.lines(x=unit(c(13.5,23.5),'mm'),y=unit(c(8,8),'mm'));grid.text('10 mm',x=unit(18.5,'mm'),y=unit(5,'mm'),gp=gpar(fontfamily='Arial',fontsize=6));grDevices::dev.off()
  grDevices::cairo_pdf(file.path(outdir,'ACTUAL_SIZE_PROOFS',paste0(fig,'_PANEL_PROOFS.pdf')),width=210/25.4,height=297/25.4,family='Arial')
  for(i in seq_len(nrow(ss))){si<-ss[i,];l<-coords[coords$panel_id==si$panel_id,];grid.newpage();grid.text(paste(si$panel_id,'| print at 100%'),x=unit(12,'mm'),y=unit(285,'mm'),just=c('left','top'),gp=gpar(fontfamily='Arial',fontsize=8));draw_panel(si$panel_id,x=12,y=20,w=as.numeric(l$width_mm),h=as.numeric(l$height_mm),page_h=297)}
  grDevices::dev.off();cat('COMPLETE',fig,'\n');flush.console()
}
writeLines(capture.output(sessionInfo()),file.path(outdir,'QC','R_SESSION_INFO.txt'));grDevices::dev.off();unlink(metric_file)
