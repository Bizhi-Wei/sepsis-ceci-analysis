# Corrected received-fluid associations and 72h CECI paths; observational analyses.
suppressPackageStartupMessages({library(tidyverse);library(survival)})
if (!nzchar(Sys.getenv("CECI_WORK_DIR"))) stop("Set CECI_WORK_DIR to a private working folder.")
base_dir <- normalizePath(Sys.getenv("CECI_WORK_DIR"), mustWork = TRUE)
proc_dir <- file.path(base_dir,"01_data/processed")
tab_dir <- file.path(base_dir,"03_results/tables")
cohort <- read_csv(file.path(proc_dir,"mimic_sepsis_tte_trajectory_cohort.csv"),show_col_types=FALSE)
cohort$ceci_high <- as.integer(cohort$ceci_z>=0)
cohort$hypoalbuminemia <- as.integer(cohort$albumin_24h<2.5)
covars<-c("age","sex","race_group","charlson_approx","septic_shock_icd","vasopressor_24h","imv_24h","creatinine_24h","map_min_24h")
for(v in covars) if(is.numeric(cohort[[v]]) && anyNA(cohort[[v]])) cohort[[v]][is.na(cohort[[v]])]<-median(cohort[[v]],na.rm=TRUE)
covar_form_str<-paste(covars,collapse="+")
cohort$predom_balanced_d5lr<-ifelse(cohort$balanced_vol_24h+cohort$saline_vol_24h+cohort$d5lr_vol_24h>0,as.integer((cohort$balanced_vol_24h+cohort$d5lr_vol_24h)/(cohort$balanced_vol_24h+cohort$saline_vol_24h+cohort$d5lr_vol_24h)>=.5),NA_integer_)
rows<-list();balances<-list()
for(exposure in c("predom_balanced","tte_fluid_treat","early_albumin","predom_balanced_d5lr")) {
 x<-cohort[!is.na(cohort[[exposure]]),];x$treatment<-as.integer(x[[exposure]])
 psmod<-glm(as.formula(paste("treatment~",covar_form_str)),data=x,family=binomial())
 ps<-predict(psmod,newdata=x,type="response")
 stopifnot(length(ps)==nrow(x),!anyNA(ps),all(ps>0 & ps<1))
 pt<-mean(x$treatment);wt<-ifelse(x$treatment==1,pt/ps,(1-pt)/(1-ps))
 x$wt<-pmin(pmax(wt,quantile(wt,.01)),quantile(wt,.99))
 fits<-list(Crude=coxph(Surv(time_post24h_27d,mort_28d)~treatment,data=x),Adjusted=coxph(as.formula(paste("Surv(time_post24h_27d,mort_28d)~treatment+",covar_form_str)),data=x),IPTW=coxph(Surv(time_post24h_27d,mort_28d)~treatment,data=x,weights=wt,robust=TRUE))
 inter<-summary(coxph(Surv(time_post24h_27d,mort_28d)~treatment*ceci_high,data=x,weights=wt,robust=TRUE))$coefficients["treatment:ceci_high","Pr(>|z|)"]
 for(label in c(names(fits),"IPTW_Low_CECI","IPTW_High_CECI")) {
  z<-if(label=="IPTW_Low_CECI")x[x$ceci_high==0,] else if(label=="IPTW_High_CECI")x[x$ceci_high==1,] else x
  fit<-if(label %in% names(fits)) fits[[label]] else coxph(Surv(time_post24h_27d,mort_28d)~treatment,data=z,weights=wt,robust=TRUE)
  ss<-summary(fit);ci<-ss$conf.int["treatment",]
  rows[[length(rows)+1]]<-data.frame(Exposure=exposure,Model=label,N=nrow(z),Deaths=sum(z$mort_28d),Control_N=sum(z$treatment==0),Treated_N=sum(z$treatment==1),Control_Mortality_Pct=100*mean(z$mort_28d[z$treatment==0]),Treated_Mortality_Pct=100*mean(z$mort_28d[z$treatment==1]),HR=ci[1],HR_Lower=ci[3],HR_Upper=ci[4],P_Value=ss$coefficients["treatment","Pr(>|z|)"],P_Interaction=if(label=="IPTW")inter else NA)
 }
 mm<-model.matrix(psmod)[,-1,drop=FALSE]
 for(v in colnames(mm)) {
  val<-mm[,v];den<-sqrt((var(val[x$treatment==1])+var(val[x$treatment==0]))/2)
  smd<-function(w)abs(weighted.mean(val[x$treatment==1],w[x$treatment==1])-weighted.mean(val[x$treatment==0],w[x$treatment==0]))/den
  balances[[length(balances)+1]]<-data.frame(Exposure=exposure,Covariate=v,Unweighted_SMD=smd(rep(1,nrow(x))),Weighted_SMD=smd(x$wt))
 }
}
write_csv(bind_rows(rows),file.path(tab_dir,"Table_S11_Exploratory_Fluid_Associations.csv"))
write_csv(bind_rows(balances),file.path(tab_dir,"Table_S11b_Fluid_Weighting_Balance.csv"))
write_csv(bind_rows(rows) %>% filter(Model=="IPTW",Exposure %in% c("predom_balanced","tte_fluid_treat","early_albumin")),file.path(tab_dir,"Table_7_Exploratory_Fluid_Summary.csv"))
# PACKAGE 2: 72-HOUR SEQUENTIAL LANDMARK & MULTI-STATE TRAJECTORIES
# ==============================================================================
cat("\n=== [3/4] Package 2: 72-Hour Sequential Landmark Dynamic Modeling ===\n")

# Sequential Landmark Cohorts
# 1. 24h survivors (N=3,254 baseline complete cases, follow-up starting at +24h)
cohort_24h <- cohort %>% filter(time_28d > 1.0)
# 2. 48h survivors (surviving > 48h from admission, i.e., time_post24h_27d > 1.0)
cohort_48h <- cohort %>% filter(time_post24h_27d > 1.0) %>% mutate(time_post48h = pmin(time_post24h_27d - 1.0, 26.0))
# 3. 72h survivors (surviving > 72h from admission, i.e., time_post24h_27d > 2.0)
cohort_72h <- cohort %>% filter(time_post24h_27d > 2.0) %>% mutate(time_post72h = pmin(time_post24h_27d - 2.0, 25.0))

cat(sprintf("24h Landmark Cohort: N = %d (Deaths = %d, %.1f%%)\n", nrow(cohort_24h), sum(cohort_24h$mort_28d), mean(cohort_24h$mort_28d) * 100))
cat(sprintf("48h Landmark Cohort: N = %d (Deaths = %d, %.1f%%)\n", nrow(cohort_48h), sum(cohort_48h$mort_28d), mean(cohort_48h$mort_28d) * 100))
cat(sprintf("72h Landmark Cohort: N = %d (Deaths = %d, %.1f%%)\n", nrow(cohort_72h), sum(cohort_72h$mort_28d), mean(cohort_72h$mort_28d) * 100))

# Fit Landmark Cox models
cox_lm24 <- coxph(Surv(time_post24h_27d, mort_28d) ~ ceci_high + age + sex + charlson_approx + septic_shock_icd + vasopressor_24h + imv_24h + creatinine_24h, data = cohort_24h)
cox_lm48 <- coxph(Surv(time_post48h, mort_28d) ~ I(ceci_d2_z >= 0) + age + sex + charlson_approx + septic_shock_icd + vasopressor_24h + imv_24h + creatinine_24h, data = cohort_48h)
cox_lm72 <- coxph(Surv(time_post72h, mort_28d) ~ I(ceci_d3_z >= 0) + age + sex + charlson_approx + septic_shock_icd + vasopressor_24h + imv_24h + creatinine_24h, data = cohort_72h)

# 8 Distinct 3-Day Multi-State Trajectories in 72h Survivors (N = 3,005)
cohort_72h <- cohort_72h %>%
  mutate(
    s1 = ifelse(ceci_z >= 0, "High", "Low"),
    s2 = ifelse(ceci_d2_z >= 0, "High", "Low"),
    s3 = ifelse(ceci_d3_z >= 0, "High", "Low"),
    trajectory_3d = factor(
      paste(s1, "->", s2, "->", s3),
      levels = c(
        "Low -> Low -> Low",
        "High -> Low -> Low",
        "High -> High -> Low",
        "Low -> High -> Low",
        "High -> Low -> High",
        "Low -> High -> High",
        "Low -> Low -> High",
        "High -> High -> High"
      )
    )
  )

cox_traj_8state <- coxph(
  Surv(time_post72h, mort_28d) ~ trajectory_3d + age + sex + charlson_approx + septic_shock_icd + vasopressor_24h + imv_24h + creatinine_24h,
  data = cohort_72h
)
sum_traj_8state <- summary(cox_traj_8state)

# Build Table S8
traj_levels <- levels(cohort_72h$trajectory_3d)
table_s8_rows <- list()

for (i in seq_along(traj_levels)) {
  lvl <- traj_levels[i]
  sub_df <- cohort_72h %>% filter(trajectory_3d == lvl)
  n_pts <- nrow(sub_df)
  d_pts <- sum(sub_df$mort_28d)
  mort_pct <- round(d_pts / n_pts * 100, 1)
  
  if (i == 1) {
    # Reference group
    table_s8_rows[[i]] <- tibble(
      Trajectory_Path = lvl,
      Trajectory_Class = "Persistent Low Injury (Reference)",
      N_Patients = n_pts,
      Deaths_Post72h = d_pts,
      Mortality_Pct = mort_pct,
      Adjusted_HR = 1.00,
      CI_95_Lower = 1.00,
      CI_95_Upper = 1.00,
      P_Value = NA_real_
    )
  } else {
    coef_name <- paste0("trajectory_3d", lvl)
    cf <- sum_traj_8state$coefficients[coef_name, ]
    ci <- sum_traj_8state$conf.int[coef_name, ]
    
    t_class <- case_when(
      lvl == "High -> High -> High" ~ "Persistently High CECI",
      lvl %in% c("High -> Low -> Low", "High -> High -> Low") ~ "Decreasing CECI",
      lvl %in% c("Low -> High -> High", "Low -> Low -> High", "Low -> High -> Low") ~ "Late or transient high CECI",
      TRUE ~ "Fluctuating CECI"
    )
    
    table_s8_rows[[i]] <- tibble(
      Trajectory_Path = lvl,
      Trajectory_Class = t_class,
      N_Patients = n_pts,
      Deaths_Post72h = d_pts,
      Mortality_Pct = mort_pct,
      Adjusted_HR = round(ci["exp(coef)"], 2),
      CI_95_Lower = round(ci["lower .95"], 2),
      CI_95_Upper = round(ci["upper .95"], 2),
      P_Value = signif(cf["Pr(>|z|)"], 3)
    )
  }
}

table_s8 <- bind_rows(table_s8_rows)

# Export Table S8
write_csv(table_s8, file.path(tab_dir, "Table_S12_72Hour_CECI_Paths.csv"))

cat("Exported: Table_S12_72Hour_CECI_Paths.csv\n")
print(table_s8)

cat("\n=== [4/4] Finished 04b_fluid_tte_and_longitudinal_landmarks.R Successfully ===\n")
