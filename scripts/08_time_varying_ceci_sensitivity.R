# Post hoc PH sensitivity: time origin is ICU Hour 24, cut = 6 days since
# landmark, corresponding to ICU Day 7. Events exactly at Day 7 belong early.
# Existing primary adjustment set is retained, with period-specific baseline
# hazards and CECI coefficients and patient-clustered robust variance.
suppressPackageStartupMessages({library(tidyverse);library(survival)})
if (!nzchar(Sys.getenv("CECI_WORK_DIR"))) stop("Set CECI_WORK_DIR to a private working folder.")
P <- normalizePath(Sys.getenv("CECI_WORK_DIR"), mustWork = TRUE)
O <- file.path(P,'03_results/tables')
a <- readr::read_csv(file.path(P,'01_data/processed/mimic_sepsis_endothelial_cohort.csv'),show_col_types=FALSE)
a <- a %>% mutate(
 bilirubin_imp=ifelse(is.na(bilirubin_24h),median(bilirubin_24h,na.rm=TRUE),bilirubin_24h),
 bilirubin_na=as.integer(is.na(bilirubin_24h)),
 wbc_imp=ifelse(is.na(wbc_24h),median(wbc_24h,na.rm=TRUE),wbc_24h),
 map_mean_imp=ifelse(is.na(map_mean_24h),median(map_mean_24h,na.rm=TRUE),map_mean_24h),
 careunit_group=case_when(grepl('Medical',first_careunit)~'MICU',grepl('Surgical|Trauma',first_careunit)~'SICU/TSICU',grepl('Coronary|Cardiac',first_careunit)~'CCU/CVICU',TRUE~'Other_ICU'),
 race_clean=case_when(race_group=='White'~'White',race_group=='Black'~'Black',TRUE~'Other'))
covars <- c('age','sex','race_clean','vasopressor_24h','imv_24h','rrt_24h','careunit_group',
 'charlson_approx','cirrhosis_portal_htn','ckd','chf','diabetes','malignancy',
 'source_pulmonary','source_abdominal','source_urinary','creatinine_24h','map_mean_imp',
 'bilirubin_imp','bilirubin_na','wbc_imp')
stopifnot(all(a$time_post24h_27d>0),all(complete.cases(a[,covars])),a$stay_id %>% n_distinct()==nrow(a))
primary<-coxph(as.formula(paste('Surv(time_post24h_27d,mort_28d)~ceci_z+',paste(covars,collapse='+'))),data=a)
ph<-cox.zph(primary)
write.csv(data.frame(term=rownames(ph$table),ph$table),file.path(O,'Table_S16_Proportional_Hazards_Diagnostics.csv'),row.names=FALSE)
long<-survSplit(Surv(time_post24h_27d,mort_28d)~.,data=a,cut=6,episode='period',id='split_patient')
long$ceci_early<-long$ceci_z*as.integer(long$period==1)
long$ceci_late<-long$ceci_z*as.integer(long$period==2)
stopifnot(all(long$time_post24h_27d[long$period==1]<=6),all(long$tstart[long$period==2]==6),
 sum(long$mort_28d)==sum(a$mort_28d),nrow(long[long$period==2,])==sum(a$time_post24h_27d>6))
form<-as.formula(paste('Surv(tstart,time_post24h_27d,mort_28d)~ceci_early+ceci_late+strata(period)+',paste(covars,collapse='+')))
fit<-coxph(form,data=long,cluster=stay_id,robust=TRUE)
sm<-summary(fit);b<-coef(fit);v<-vcov(fit)
contrast<-b['ceci_early']-b['ceci_late']
contrast_se<-sqrt(v['ceci_early','ceci_early']+v['ceci_late','ceci_late']-2*v['ceci_early','ceci_late'])
interaction_p<-2*pnorm(abs(contrast/contrast_se),lower.tail=FALSE)
rows<-lapply(1:2,function(k){
 term<-c('ceci_early','ceci_late')[k];r<-long[long$period==k,]
 data.frame(interval=c('ICU Hour 24 through Day 7 inclusive','After ICU Day 7 through Day 28')[k],
  start_days_since_24h=c(0,6)[k],end_days_since_24h=c(6,27)[k],
  people_at_interval_start=nrow(r),events=sum(r$mort_28d),person_days=sum(r$time_post24h_27d-r$tstart),
  HR=exp(b[term]),CI95_low=exp(b[term]-1.96*sqrt(v[term,term])),CI95_high=exp(b[term]+1.96*sqrt(v[term,term])),
  robust_SE=sqrt(v[term,term]),P=2*pnorm(abs(b[term]/sqrt(v[term,term])),lower.tail=FALSE),
  period_CECI_interaction_Wald_P=interaction_p,
  method='Original full adjustment; period strata and CECI coefficients; patient-clustered robust variance; post hoc sensitivity')
})
result<-bind_rows(rows)
write.csv(result,file.path(O,'Table_S17_CECI_Time_Varying_Associations.csv'),row.names=FALSE)
write.csv(data.frame(metric=c('unique_patients','original_events','split_records','events_exactly_ICU_day7','events_early','events_late','patients_at_hour24','patients_at_day7','CECI_Schoenfeld_P','global_Schoenfeld_P','CECI_period_interaction_Wald_P'),value=c(nrow(a),sum(a$mort_28d),nrow(long),sum(a$mort_28d==1 & a$time_28d==7),sum(long$mort_28d[long$period==1]),sum(long$mort_28d[long$period==2]),nrow(a),sum(a$time_28d>7),ph$table['ceci_z','p'],ph$table['GLOBAL','p'],interaction_p)),file.path(O,'Table_S17b_Time_Varying_Risksets.csv'),row.names=FALSE)
print(result);print(ph$table[c('ceci_z','GLOBAL'),,drop=FALSE])
