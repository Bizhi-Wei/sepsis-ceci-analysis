options(warn=1)
suppressPackageStartupMessages({library(pROC);library(ggplot2);library(cowplot);library(readr)})
if (!nzchar(Sys.getenv("CECI_WORK_DIR"))) stop("Set CECI_WORK_DIR.")
O <- file.path(normalizePath(Sys.getenv("CECI_WORK_DIR"), mustWork = TRUE), "eicu_analysis")
x<-read.csv(file.path(O,'restricted_local_only/eicu_complete_frozen_ceci.csv'))
gate<-jsonlite::fromJSON(file.path(O,'extraction_verification.json'))
stopifnot(gate$score_association_assessment_possible,nrow(x)==4169,!anyNA(x[,c('ceci_frozen','age','sex','hospital_mortality')]),!anyDuplicated(x$uniquepid))
y<-x$hospital_mortality;s<-x$ceci_frozen
hospital<-sort(unique(x$hospitalid));site<-match(x$hospitalid,hospital)
ng<-length(hospital)
ord<-order(s);ys<-y[ord];ss<-s[ord]
tie<-cumsum(c(TRUE,diff(ss)!=0))
aucw<-function(w){
 w<-w[ord];sums<-rowsum(cbind(w*ys,w*(1-ys)),tie,reorder=FALSE)
 np<-sum(sums[,1]);nn<-sum(sums[,2]);if(np==0 || nn==0)return(NA_real_)
 sum(sums[,1]*(cumsum(sums[,2])-sums[,2]/2))/(np*nn)
}
X0<-cbind(Intercept=1,CECI=s)
X1<-cbind(Intercept=1,CECI=s,Age=x$age,Sex=x$sex)
fitbeta<-function(X,w){
 keep<-w>0
 fit<-suppressWarnings(glm.fit(X[keep,,drop=FALSE],y[keep],weights=w[keep],family=binomial()))
 if(!fit$converged || !all(is.finite(fit$coefficients)))return(NA_real_)
 unname(fit$coefficients['CECI'])
}
metrics<-function(w){
 high<-s>=0
 c(AUC=aucw(w),logOR=fitbeta(X0,w),logOR_age_sex=fitbeta(X1,w),mortality_low=sum(w[!high]*y[!high])/sum(w[!high]),mortality_high=sum(w[high]*y[high])/sum(w[high]))
}
point<-metrics(rep(1,nrow(x)))
stopifnot(abs(point['AUC']-as.numeric(auc(roc(y,s,direction='<',quiet=TRUE))))<1e-12)
RNGkind('Mersenne-Twister','Inversion','Rejection');set.seed(20261105);B<-1000
boot<-matrix(NA_real_,nrow=B,ncol=length(point),dimnames=list(NULL,names(point)))
for(b in 1:B){
 counts<-tabulate(sample.int(ng,ng,replace=TRUE),nbins=ng)
 boot[b,]<-metrics(counts[site])
}
write.csv(data.frame(replicate=1:B,boot),file.path(O,'hospital_cluster_bootstrap_metrics.csv'),row.names=FALSE)
ci<-apply(boot,2,quantile,probs=c(.025,.975),na.rm=TRUE)
summary<-data.frame(Metric=c('Fixed CECI AUC','Unadjusted OR per MIMIC CECI SD','Age/sex adjusted OR per MIMIC CECI SD','Observed mortality: CECI below 0','Observed mortality: CECI at least 0'),Estimate=c(point[1],exp(point[2:3]),point[4:5]),CI95_Lower=c(ci[1,1],exp(ci[1,2:3]),ci[1,4:5]),CI95_Upper=c(ci[2,1],exp(ci[2,2:3]),ci[2,4:5]),Valid_Bootstrap_Replicates=colSums(is.finite(boot)),External_Endpoint='In-hospital mortality conditional on survival/stay through ICU Hour 24',Method='1000 hospital-cluster bootstrap resamples; frozen MIMIC laboratory transformation; seed 20261105')
write.csv(summary,file.path(O,'external_score_assessment.csv'),row.names=FALSE)
groups<-do.call(rbind,lapply(c(FALSE,TRUE),function(h){z<-x[(x$ceci_frozen>=0)==h,];j<-if(h)5 else 4;data.frame(Group=if(h)'CECI >= 0' else 'CECI < 0',N=nrow(z),Hospital_deaths=sum(z$hospital_mortality),Mortality=mean(z$hospital_mortality),CI95_Lower=ci[1,j],CI95_Upper=ci[2,j])}))
write.csv(groups,file.path(O,'predefined_score_group_mortality.csv'),row.names=FALSE)
r<-roc(y,s,direction='<',quiet=TRUE)
curve<-data.frame(FPR=1-r$specificities,TPR=r$sensitivities)
write.csv(curve,file.path(O,'external_score_ROC_source.csv'),row.names=FALSE)
theme_set(theme_classic(base_size=10,base_family='sans'))
a<-ggplot(curve,aes(FPR,TPR))+geom_abline(slope=1,intercept=0,lty=2,color='grey70')+geom_path(linewidth=.85,color='#2166AC')+coord_equal()+labs(title='A  Fixed CECI discrimination',subtitle=sprintf('AUC %.3f (95%% CI %.3f to %.3f)',point[1],ci[1,1],ci[2,1]),x='False-positive rate',y='True-positive rate')+theme(plot.title=element_text(face='bold',size=11),plot.subtitle=element_text(size=9))
groups$Group<-factor(groups$Group,levels=c('CECI < 0','CECI >= 0'))
b<-ggplot(groups,aes(Group,100*Mortality))+geom_col(fill='#2166AC',width=.55)+geom_errorbar(aes(ymin=100*CI95_Lower,ymax=100*CI95_Upper),width=.12,linewidth=.5)+geom_text(aes(y=100*CI95_Upper+2,label=sprintf('%d/%d deaths',Hospital_deaths,N)),size=3)+scale_y_continuous(limits=c(0,max(100*groups$CI95_Upper)+8),expand=expansion(mult=c(0,.03)))+labs(title='B  Observed hospital mortality',subtitle='MIMIC scaling and zero cutoff; observed group rates',x=NULL,y='Hospital mortality (%)')+theme(plot.title=element_text(face='bold',size=11),plot.subtitle=element_text(size=9))
fig<-plot_grid(a,b,ncol=2,rel_widths=c(1,1.05))
ggsave(file.path(O,'External_CECI_Hospital_Outcome_Assessment.pdf'),fig,width=10,height=4.5,device=cairo_pdf,bg='white')
ggsave(file.path(O,'External_CECI_Hospital_Outcome_Assessment.png'),fig,width=10,height=4.5,dpi=300,bg='white')
jsonlite::write_json(list(N=nrow(x),events=sum(y),hospitals=ng,replicates=B,valid_replicates=colSums(is.finite(boot)),point_metrics=point,reference='Outcome is hospital mortality, not 28-day vital status. OR models refit in eICU are association models, not locked augmented-model validation.',ceci_frozen_mean=mean(s),ceci_frozen_sd=sd(s),strict_culture_n=sum(x$culture_recorded),strict_culture_events=sum(y[x$culture_recorded==1]),strict_culture_inference='Descriptive only: 88 people/19 events, sparse interface and selection'),file.path(O,'assessment_verification.json'),pretty=TRUE,auto_unbox=TRUE)
print(summary);print(groups)
