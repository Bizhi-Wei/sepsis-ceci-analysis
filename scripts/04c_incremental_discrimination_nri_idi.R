# Corrected internal clinical validation and observational sensitivity analyses.
suppressPackageStartupMessages({library(survival);library(pROC);library(readr)})
if (!nzchar(Sys.getenv("CECI_WORK_DIR"))) stop("Set CECI_WORK_DIR to a private working folder.")
P<-normalizePath(Sys.getenv("CECI_WORK_DIR"), mustWork = TRUE)
tab_dir<-file.path(P,"03_results/tables")
a<-read.csv(file.path(P,"01_data/processed/mimic_sepsis_tte_trajectory_cohort.csv"))
# Internal ten-fold predictions; every model and score-scaling parameter is fit in training only.
set.seed(20261104)
fold<-integer(nrow(a))
for(y in 0:1){ix<-which(a$mort_28d==y);fold[ix]<-sample(rep(1:10,length.out=length(ix)))}
basevars<-c('age','sex','charlson_approx','vasopressor_24h','imv_24h','rrt_24h')
pred0<-pred1<-numeric(nrow(a))
rawscore<-function(x) cbind(log2(x$lactate_24h),-x$albumin_24h,-log2(x$platelet_24h),log2(x$inr_24h))
for(k in 1:10) {
 train<-a[fold!=k,];test<-a[fold==k,]
 raw<-rawscore(train);mu<-colMeans(raw);sdv<-apply(raw,2,sd)
 score<-function(x)rowMeans(sweep(sweep(rawscore(x),2,mu,'-'),2,sdv,'/'))
 train$ceci_cv<-score(train);test$ceci_cv<-score(test)
 b<-glm(as.formula(paste('mort_28d~',paste(basevars,collapse='+'))),data=train,family=binomial())
 aug<-update(b,.~.+ceci_cv)
 pred0[fold==k]<-predict(b,newdata=test,type='response');pred1[fold==k]<-predict(aug,newdata=test,type='response')
}
out<-list()
for(nm in c('Clinical base','Clinical base plus CECI')) {
 pr<-if(nm=='Clinical base')pred0 else pred1
 rocobj<-roc(a$mort_28d,pr,direction='<',quiet=TRUE)
 lp<-qlogis(pr);cal<-glm(a$mort_28d~lp,family=binomial())
 out[[length(out)+1]]<-data.frame(model=nm,validation='Stratified ten-fold out-of-fold; seed 20261104; CECI scales fit in each training fold',n=nrow(a),events=sum(a$mort_28d),AUC=as.numeric(auc(rocobj)),Brier=mean((pr-a$mort_28d)^2),calibration_intercept=coef(cal)[1],calibration_slope=coef(cal)[2])
}
write.csv(do.call(rbind,out),file.path(tab_dir,'Table_S13_Clinical_Internal_Validation.csv'),row.names=FALSE)
saveRDS(data.frame(outcome=a$mort_28d,fold=fold,clinical_base_probability=pred0,ceci_augmented_probability=pred1),file.path(P,'01_data/processed/clinical_internal_oof_predictions.rds'))

# Common-outcome HR-to-RR approximation from VanderWeele and Ding (2017).
iptw<-read.csv(file.path(P,'03_results/tables/Table_3_MIMIC_IPTW_Association.csv'))
tr48<-read.csv(file.path(P,'03_results/tables/Table_6_Dynamic_CECI_Trajectory_Outcomes.csv'))
tr72<-read.csv(file.path(P,'03_results/tables/Table_S12_72Hour_CECI_Paths.csv'))
vals<-data.frame(analysis=c('High versus low CECI','48h persistent high versus persistent low','72h High-High-High versus Low-Low-Low'),HR=c(iptw$HR,tr48$Adjusted_HR[4],tr72$Adjusted_HR[8]),lower=c(iptw$HR_Lower,tr48$CI_95_Low[4],tr72$CI_95_Lower[8]))
hr2rr<-function(h)(1-.5^sqrt(h))/(1-.5^sqrt(1/h))
ev<-function(r)ifelse(r>1,r+sqrt(r*(r-1)),1)
vals$RR_approx<-hr2rr(vals$HR);vals$Evalue_common_outcome<-ev(vals$RR_approx);vals$Evalue_lower_common_outcome<-ev(hr2rr(vals$lower));vals$Evalue_direct_HR<-ev(vals$HR)
write.csv(vals,file.path(tab_dir,'Table_S14_Evalue_Common_Outcome_Sensitivity.csv'),row.names=FALSE)

suppressPackageStartupMessages(library(dplyr))
df_full<-readr::read_csv(file.path(P,"01_data/processed/mimic_sepsis_eligible_full.csv"),show_col_types=FALSE)
# 2. Missingness IPCW / IPMW across Full Eligible Cohort (N = 7,195)
# Fit multivariable logistic regression for being a complete case
m_missing <- glm(
  complete_panel_24h ~ age + sex + septic_shock_icd + vasopressor_24h + imv_24h +
    cirrhosis_portal_htn + chf + ckd + diabetes + malignancy,
  data = df_full,
  family = binomial()
)

df_full$p_obs <- predict(m_missing, newdata = df_full, type = "response")
p_marginal <- mean(df_full$complete_panel_24h)
df_full$ipmw <- ifelse(
  df_full$complete_panel_24h == 1,
  p_marginal / df_full$p_obs,
  (1 - p_marginal) / (1 - df_full$p_obs)
)

# Extract weights for complete cases (analytic cohort)
complete_weights <- df_full %>%
  filter(complete_panel_24h == 1) %>%
  select(stay_id, ipmw)

cohort_ipcw <- a %>%
  left_join(complete_weights, by = "stay_id")

wt_mean <- mean(cohort_ipcw$ipmw)
wt_sd   <- sd(cohort_ipcw$ipmw)
wt_min  <- min(cohort_ipcw$ipmw)
wt_p01  <- quantile(cohort_ipcw$ipmw, 0.01)
wt_p50  <- median(cohort_ipcw$ipmw)
wt_p99  <- quantile(cohort_ipcw$ipmw, 0.99)
wt_max  <- max(cohort_ipcw$ipmw)

cat(sprintf("Stabilized IPMW Weights: Mean = %.4f (SD = %.3f), Median = %.3f, Range = [%.3f, %.3f], 99th Pct = %.3f\n",
            wt_mean, wt_sd, wt_p50, wt_min, wt_max, wt_p99))

# Multivariable Cox: Complete Case Unweighted vs IPCW Weighted
cov_formula <- "age + sex + race_group + charlson_approx + septic_shock_icd + vasopressor_24h + imv_24h + creatinine_24h"

cox_unweighted <- coxph(as.formula(paste("Surv(time_28d, mort_28d) ~ ceci_z +", cov_formula)), data = cohort_ipcw)
cox_ipcw       <- coxph(as.formula(paste("Surv(time_28d, mort_28d) ~ ceci_z +", cov_formula)),
                        data = cohort_ipcw, weights = ipmw, robust = TRUE)

sum_unwt <- summary(cox_unweighted)
sum_ipcw <- summary(cox_ipcw)

hr_unwt <- sum_unwt$coefficients["ceci_z", "exp(coef)"]
ci_unwt_lo <- sum_unwt$conf.int["ceci_z", "lower .95"]
ci_unwt_hi <- sum_unwt$conf.int["ceci_z", "upper .95"]
p_unwt <- sum_unwt$coefficients["ceci_z", "Pr(>|z|)"]

hr_ipcw <- sum_ipcw$coefficients["ceci_z", "exp(coef)"]
ci_ipcw_lo <- sum_ipcw$conf.int["ceci_z", "lower .95"]
ci_ipcw_hi <- sum_ipcw$conf.int["ceci_z", "upper .95"]
p_ipcw <- sum_ipcw$coefficients["ceci_z", "Pr(>|z|)"]

pct_diff_log_hr <- abs(log(hr_ipcw) - log(hr_unwt)) / log(hr_unwt) * 100


write.csv(data.frame(Model=c("Complete cases","IPMW weighted complete cases"),N=nrow(cohort_ipcw),HR=c(hr_unwt,hr_ipcw),HR_Lower=c(ci_unwt_lo,ci_ipcw_lo),HR_Upper=c(ci_unwt_hi,ci_ipcw_hi),P_Value=c(p_unwt,p_ipcw),Selection_Population_N=nrow(df_full),Mean_Weight=c(1,wt_mean),P99_Weight=c(1,wt_p99)),file.path(tab_dir,"Table_S15_Missingness_Weighting_Sensitivity.csv"),row.names=FALSE)
