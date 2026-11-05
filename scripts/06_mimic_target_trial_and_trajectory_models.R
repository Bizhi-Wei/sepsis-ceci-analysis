# ==============================================================================
# 06_mimic_target_trial_and_trajectory_models.R
# Target Trial Emulation (TTE) of Fluid Resuscitation Intensity & 
# 72-Hour Dynamic Endothelial Injury Trajectories in MIMIC-IV Sepsis Cohort
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(survival)
  library(ggplot2)
  library(cowplot)
  library(gridExtra)
})

if (!nzchar(Sys.getenv("CECI_WORK_DIR"))) stop("Set CECI_WORK_DIR to a private working folder.")
base_dir <- normalizePath(Sys.getenv("CECI_WORK_DIR"), mustWork = TRUE)
proc_dir <- file.path(base_dir, "01_data", "processed")
tab_dir <- file.path(base_dir, "03_results", "tables")
fig_dir <- file.path(base_dir, "03_results", "figures")

dir.create(tab_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

cat("=== 1. Loading Augmented TTE & Trajectory Sepsis Cohort ===\n")
cohort <- read_csv(file.path(proc_dir, "mimic_sepsis_tte_trajectory_cohort.csv"), show_col_types = FALSE)
cat(sprintf("Cohort size: %d patients\n", nrow(cohort)))

# Define baseline high CECI indicator (standard zero z-score threshold matching Table 1-3)
cohort <- cohort %>%
  mutate(
    ceci_high = ifelse(ceci_z >= 0, 1, 0),
    ceci_endotype = ifelse(ceci_z >= 0, "High_Endothelial_Injury", "Low_Endothelial_Injury"),
    ceci_trajectory_group = factor(ceci_trajectory_group, 
      levels = c("Persistent_Low_Injury", "Resolving_Injury", "Deteriorating_Injury", "Persistent_High_Injury")),
    tte_fluid_group = factor(tte_fluid_group, 
      levels = c("Restrictive_lt_30ml_kg", "Liberal_ge_30ml_kg"))
  )

# ==============================================================================
# 2. 48-Hour Landmark Dynamic Trajectory Modeling
# ==============================================================================
cat("\n=== 2. 48-Hour Landmark Dynamic Trajectory Survival Analysis ===\n")
# To strictly eliminate immortal time bias, trajectory groups are evaluated
# among patients surviving past the completion of Day 2 (time_post24h_27d > 1.0 day)
cohort_48h <- cohort %>%
  filter(time_post24h_27d > 1.0) %>%
  mutate(
    time_post48h_26d = pmin(time_post24h_27d - 1.0, 26.0),
    ceci_trajectory_group = factor(ceci_trajectory_group, 
      levels = c("Persistent_Low_Injury", "Resolving_Injury", "Deteriorating_Injury", "Persistent_High_Injury"))
  )
cat(sprintf("48h Landmark Cohort size: %d patients (%d excluded early deaths between 24h and 48h)\n", 
            nrow(cohort_48h), nrow(cohort) - nrow(cohort_48h)))

# Kaplan-Meier survival analysis across dynamic trajectories in 48h landmark cohort
fit_traj <- survfit(Surv(time_post48h_26d, mort_28d) ~ ceci_trajectory_group, data = cohort_48h)
surv_diff <- survdiff(Surv(time_post48h_26d, mort_28d) ~ ceci_trajectory_group, data = cohort_48h)
logrank_p <- pchisq(surv_diff$chisq, length(surv_diff$n) - 1, lower.tail = FALSE)
cat(sprintf("Log-rank test P-value across 4 trajectories: %.2e\n", logrank_p))

# Multivariable Cox regression for dynamic trajectory groups
cox_traj <- coxph(
  Surv(time_post48h_26d, mort_28d) ~ ceci_trajectory_group + age + sex + charlson_approx + 
    septic_shock_icd + vasopressor_24h + imv_24h + creatinine_24h,
  data = cohort_48h
)
cox_traj_sum <- summary(cox_traj)

traj_table <- tibble(
  Trajectory_Group = levels(cohort_48h$ceci_trajectory_group),
  N_Patients = as.numeric(table(cohort_48h$ceci_trajectory_group)),
  Deaths_28d = as.numeric(tapply(cohort_48h$mort_28d, cohort_48h$ceci_trajectory_group, sum)),
  Mortality_Rate_Pct = round(Deaths_28d / N_Patients * 100, 1),
  Adjusted_HR = c(1.0, 
                  round(cox_traj_sum$coefficients[1, "exp(coef)"], 2),
                  round(cox_traj_sum$coefficients[2, "exp(coef)"], 2),
                  round(cox_traj_sum$coefficients[3, "exp(coef)"], 2)),
  CI_95_Low = c(1.0, 
                round(cox_traj_sum$conf.int[1, "lower .95"], 2),
                round(cox_traj_sum$conf.int[2, "lower .95"], 2),
                round(cox_traj_sum$conf.int[3, "lower .95"], 2)),
  CI_95_High = c(1.0, 
                 round(cox_traj_sum$conf.int[1, "upper .95"], 2),
                 round(cox_traj_sum$conf.int[2, "upper .95"], 2),
                 round(cox_traj_sum$conf.int[3, "upper .95"], 2)),
  P_Value = c(NA, 
              signif(cox_traj_sum$coefficients[1, "Pr(>|z|)"], 3),
              signif(cox_traj_sum$coefficients[2, "Pr(>|z|)"], 3),
              signif(cox_traj_sum$coefficients[3, "Pr(>|z|)"], 3))
)

write_csv(traj_table, file.path(tab_dir, "Table_6_Dynamic_CECI_Trajectory_Outcomes.csv"))
cat("Saved: Table_6_Dynamic_CECI_Trajectory_Outcomes.csv\n")
print(traj_table)

if ('--tables-only' %in% commandArgs(trailingOnly=TRUE)) quit(status=0)
# Panel A: 72-Hour Longitudinal CECI Trajectory Plot
traj_counts <- table(cohort_48h$ceci_trajectory_group)
label_low <- sprintf("Persistent Low CECI (N=%d)", traj_counts["Persistent_Low_Injury"])
label_res <- sprintf("Resolving CECI (N=%d)", traj_counts["Resolving_Injury"])
label_det <- sprintf("Deteriorating CECI (N=%d)", traj_counts["Deteriorating_Injury"])
label_high <- sprintf("Persistent High CECI (N=%d)", traj_counts["Persistent_High_Injury"])

traj_summary <- cohort_48h %>%
  group_by(ceci_trajectory_group) %>%
  summarise(
    D1_mean = mean(ceci_z, na.rm = TRUE),
    D1_se = sd(ceci_z, na.rm = TRUE) / sqrt(n()),
    D2_mean = mean(ceci_d2_z, na.rm = TRUE),
    D2_se = sd(ceci_d2_z, na.rm = TRUE) / sqrt(n()),
    D3_mean = mean(ceci_d3_z, na.rm = TRUE),
    D3_se = sd(ceci_d3_z, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  ) %>%
  pivot_longer(cols = starts_with("D"), 
               names_to = c("Day", ".value"), 
               names_pattern = "(D[1-3])_(mean|se)") %>%
  mutate(
    Day_Num = case_when(Day == "D1" ~ 1, Day == "D2" ~ 2, Day == "D3" ~ 3),
    Trajectory = case_when(
      ceci_trajectory_group == "Persistent_Low_Injury" ~ label_low,
      ceci_trajectory_group == "Resolving_Injury" ~ label_res,
      ceci_trajectory_group == "Deteriorating_Injury" ~ label_det,
      ceci_trajectory_group == "Persistent_High_Injury" ~ label_high
    )
  )

palette_traj <- c("#3498DB", "#2ECC71", "#E67E22", "#C0392B")
names(palette_traj) <- c(label_low, label_res, label_det, label_high)

panel_a <- ggplot(traj_summary, aes(x = Day_Num, y = mean, color = Trajectory, group = Trajectory)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = mean - 1.96 * se, ymax = mean + 1.96 * se), width = 0.08, linewidth = 0.8) +
  scale_color_manual(values = palette_traj) +
  scale_x_continuous(breaks = c(1, 2, 3), labels = c("Baseline (-6 to +24h)", "Day 2 (24-48h)", "Day 3 (48-72h)")) +
  theme_classic(base_size = 11) +
  labs(
    title = "A  CECI across fixed daily windows",
    subtitle = sprintf("Serial CECI evolution in 48-hour landmark survivors (N = %s)", format(nrow(cohort_48h), big.mark = ",")),
    x = "ICU Monitoring Period",
    y = "CECI Score (Z-Score SD)",
    color = "Dynamic Trajectory"
  ) +
  theme(plot.title = element_text(face = "bold", size = 12),
        legend.position = "right")

# Panel B: Kaplan-Meier Survival by Trajectory (48-hour Landmark Cohort)
surv_obj <- survfit(Surv(time_post48h_26d, mort_28d) ~ ceci_trajectory_group, data = cohort_48h)
surv_summary <- summary(surv_obj)

km_df <- tibble(
  time = surv_summary$time,
  surv = surv_summary$surv,
  lower = surv_summary$lower,
  upper = surv_summary$upper,
  strata = factor(gsub("ceci_trajectory_group=", "", as.character(surv_summary$strata)),
                  levels = c("Persistent_Low_Injury", "Resolving_Injury", "Deteriorating_Injury", "Persistent_High_Injury"))
)

# Add origin point (time = 0, surv = 1)
km_origin <- tibble(
  time = 0,
  surv = 1.0,
  lower = 1.0,
  upper = 1.0,
  strata = factor(c("Persistent_Low_Injury", "Resolving_Injury", "Deteriorating_Injury", "Persistent_High_Injury"),
                  levels = c("Persistent_Low_Injury", "Resolving_Injury", "Deteriorating_Injury", "Persistent_High_Injury"))
)
km_df_full <- bind_rows(km_origin, km_df) %>% arrange(strata, time)

hr_high_val <- round(cox_traj_sum$coefficients["ceci_trajectory_groupPersistent_High_Injury", "exp(coef)"], 2)
hr_high_low <- round(cox_traj_sum$conf.int["ceci_trajectory_groupPersistent_High_Injury", "lower .95"], 2)
hr_high_high <- round(cox_traj_sum$conf.int["ceci_trajectory_groupPersistent_High_Injury", "upper .95"], 2)

panel_b <- ggplot(km_df_full, aes(x = time, y = surv, color = strata)) +
  geom_step(linewidth = 1.1) +
  scale_color_manual(values = c(
    "Persistent_Low_Injury" = "#3498DB",
    "Resolving_Injury" = "#2ECC71",
    "Deteriorating_Injury" = "#E67E22",
    "Persistent_High_Injury" = "#C0392B"
  ), labels = c("Persistent Low", "Resolving", "Deteriorating", "Persistent High")) +
  theme_classic(base_size = 11) +
  labs(
    title = "B  Survival after the 48-hour landmark",
    subtitle = sprintf("48h Landmark (N=%s) | Log-rank P = %.2e | High vs Low: aHR %.2f (%.2f-%.2f)",
                       format(nrow(cohort_48h), big.mark = ","), logrank_p, hr_high_val, hr_high_low, hr_high_high),
    x = "Days from 48-Hour Landmark (Follow-up to Day 28)",
    y = "Overall Survival Probability",
    color = "Dynamic Trajectory"
  ) +
  theme(plot.title = element_text(face = "bold", size = 12),
        plot.subtitle = element_text(size = 9.5, color = "grey30"),
        legend.position = "right")


# Panel C: Fixed 72h score paths; not a transition-intensity model.
t72<-read_csv(file.path(tab_dir,"Table_S12_72Hour_CECI_Paths.csv"),show_col_types=FALSE)
t72<-t72[-1,]
t72$Trajectory_Path<-factor(t72$Trajectory_Path,levels=rev(t72$Trajectory_Path))
panel_c<-ggplot(t72,aes(x=Adjusted_HR,y=Trajectory_Path))+
 geom_vline(xintercept=1,lty=2,color="grey55")+
 geom_errorbar(aes(xmin=CI_95_Lower,xmax=CI_95_Upper),orientation="y",width=.2)+
 geom_point(size=2.5,color="#2166AC")+theme_classic(base_size=11)+
 labs(title="C  CECI paths among 72-hour survivors",subtitle="Adjusted HR and 95% CI; reference Low > Low > Low",x="Hazard ratio",y=NULL)+
 theme(plot.title=element_text(face="bold",size=12),plot.subtitle=element_text(size=9))
# Panel D: Held-out predictions, with all parameters estimated in training folds.
pr<-readRDS(file.path(proc_dir,"clinical_internal_oof_predictions.rds"))
curves<-list()
for(j in c("clinical_base_probability","ceci_augmented_probability")) {
 r<-pROC::roc(pr$outcome,pr[[j]],direction="<",quiet=TRUE)
 nm<-if(j=="clinical_base_probability") "Clinical base" else "Clinical base + CECI"
 curves[[j]]<-tibble(fpr=1-r$specificities,tpr=r$sensitivities,model=sprintf("%s (AUC %.3f)",nm,as.numeric(pROC::auc(r))))
}
panel_d<-ggplot(bind_rows(curves),aes(fpr,tpr,color=model))+
 geom_abline(slope=1,intercept=0,lty=2,color="grey65")+geom_path(linewidth=.9)+
 scale_color_manual(values=c("#777777","#2166AC"))+theme_classic(base_size=11)+
 labs(title="D  Ten-fold internal clinical prediction",subtitle="3,254 patients; 1,113 deaths; held-out predictions",x="False-positive rate",y="True-positive rate",color=NULL)+
 theme(plot.title=element_text(face="bold",size=12),plot.subtitle=element_text(size=9),legend.position=c(.66,.17),legend.text=element_text(size=9))
fig6<-plot_grid(panel_a,panel_b,panel_c,panel_d,ncol=2,align="hv")
write_csv(traj_summary,file.path(tab_dir,"Figure6A_Source_Group_Means.csv"))
write_csv(km_df_full %>% select(time,surv,lower,upper,strata),file.path(tab_dir,"Figure6B_Source_KM_Estimates.csv"))
write_csv(bind_rows(curves),file.path(tab_dir,"Figure6D_Source_ROC_Curves.csv"))
stem<-file.path(fig_dir,"Fig6_CECI_Trajectories_and_Internal_Validation")
ggsave(paste0(stem,".pdf"),fig6,width=15,height=11,device=cairo_pdf)
ggsave(paste0(stem,".png"),fig6,width=15,height=11,dpi=300)
saveRDS(list(cohort=cohort,traj_table=traj_table,cox_traj=cox_traj),file.path(proc_dir,"mimic_tte_trajectory_results_bundle.rds"))
