# ==============================================================================
# 04_mimic_cox_iptw_and_figures.R
# Real-World Clinical Epidemiology in MIMIC-IV v3.1:
# Multivariable Cox Regression, Restricted Cubic Splines, Propensity Score IPTW,
# Sensitivity Landmark Audits, and Publication Figures (Fig 3 - Fig 5)
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(survival)
  library(pROC)
  library(splines)
  library(patchwork)
})

if (!nzchar(Sys.getenv("CECI_WORK_DIR"))) stop("Set CECI_WORK_DIR to a private working folder.")
BASE_DIR <- normalizePath(Sys.getenv("CECI_WORK_DIR"), mustWork = TRUE)
DATA_DIR <- file.path(BASE_DIR, "01_data", "processed")
TAB_DIR  <- file.path(BASE_DIR, "03_results", "tables")
FIG_DIR  <- file.path(BASE_DIR, "03_results", "figures")

dir.create(TAB_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

cat("=== Step 1: Loading Real MIMIC-IV Cohorts ===\n")
df_primary <- readr::read_csv(file.path(DATA_DIR, "mimic_sepsis_endothelial_cohort.csv"), show_col_types = FALSE)
df_full    <- readr::read_csv(file.path(DATA_DIR, "mimic_sepsis_eligible_full.csv"), show_col_types = FALSE)

cat(sprintf("Primary complete-case panel: N = %d (%d deaths, %.2f%%)\n",
            nrow(df_primary), sum(df_primary$mort_28d), 100 * mean(df_primary$mort_28d)))
cat(sprintf("Full eligible sepsis cohort: N = %d (%d deaths, %.2f%%)\n",
            nrow(df_full), sum(df_full$mort_28d), 100 * mean(df_full$mort_28d)))

# Standardize variable codings
df_primary <- df_primary %>%
  dplyr::mutate(
    ceci_high = ifelse(ceci_z >= 0, 1L, 0L),
    ceci_quartile = factor(dplyr::ntile(ceci_z, 4), levels = 1:4, labels = c("Q1 (Lowest)", "Q2", "Q3", "Q4 (Highest)")),
    lar_high = ifelse(lar_24h >= stats::median(lar_24h, na.rm = TRUE), 1L, 0L),
    age_scale = age / 10,
    shock_flag = ifelse(vasopressor_24h == 1 | septic_shock_icd == 1, 1L, 0L),
    # Impute missing bilirubin, wbc, and map_mean with median
    bilirubin_imp = ifelse(is.na(bilirubin_24h), stats::median(bilirubin_24h, na.rm = TRUE), bilirubin_24h),
    bilirubin_na = ifelse(is.na(bilirubin_24h), 1L, 0L),
    wbc_imp = ifelse(is.na(wbc_24h), stats::median(wbc_24h, na.rm = TRUE), wbc_24h),
    map_mean_imp = ifelse(is.na(map_mean_24h), stats::median(map_mean_24h, na.rm = TRUE), map_mean_24h),
    # Consolidate care unit
    careunit_group = case_when(
      grepl("Medical", first_careunit) ~ "MICU",
      grepl("Surgical|Trauma", first_careunit) ~ "SICU/TSICU",
      grepl("Coronary|Cardiac", first_careunit) ~ "CCU/CVICU",
      TRUE ~ "Other_ICU"
    ),
    race_clean = case_when(
      race_group == "White" ~ "White",
      race_group == "Black" ~ "Black",
      TRUE ~ "Other"
    )
  )

# ==============================================================================
# Step 2: Table 1 - Baseline Characteristics Stratified by CECI High vs Low
# ==============================================================================
cat("\n=== Step 2: Generating Table 1 Baseline Characteristics ===\n")

make_table1_row <- function(df, var_name, is_cat = FALSE) {
  g0 <- df %>% dplyr::filter(ceci_high == 0) %>% dplyr::pull(!!sym(var_name))
  g1 <- df %>% dplyr::filter(ceci_high == 1) %>% dplyr::pull(!!sym(var_name))
  gt <- df %>% dplyr::pull(!!sym(var_name))

  if (is_cat) {
    n0 <- sum(g0 == 1, na.rm = TRUE); p0 <- 100 * mean(g0 == 1, na.rm = TRUE)
    n1 <- sum(g1 == 1, na.rm = TRUE); p1 <- 100 * mean(g1 == 1, na.rm = TRUE)
    nt <- sum(gt == 1, na.rm = TRUE); pt <- 100 * mean(gt == 1, na.rm = TRUE)
    pval <- stats::chisq.test(table(df$ceci_high, df[[var_name]]))$p.value
    tibble(
      Variable = var_name,
      Overall_N_pct_or_Median_IQR = sprintf("%d (%.1f%%)", nt, pt),
      CECI_Low_N_pct_or_Median_IQR = sprintf("%d (%.1f%%)", n0, p0),
      CECI_High_N_pct_or_Median_IQR = sprintf("%d (%.1f%%)", n1, p1),
      P_Value = pval
    )
  } else {
    m0 <- stats::median(g0, na.rm = TRUE); q0_lo <- stats::quantile(g0, 0.25, na.rm = TRUE); q0_hi <- stats::quantile(g0, 0.75, na.rm = TRUE)
    m1 <- stats::median(g1, na.rm = TRUE); q1_lo <- stats::quantile(g1, 0.25, na.rm = TRUE); q1_hi <- stats::quantile(g1, 0.75, na.rm = TRUE)
    mt <- stats::median(gt, na.rm = TRUE); qt_lo <- stats::quantile(gt, 0.25, na.rm = TRUE); qt_hi <- stats::quantile(gt, 0.75, na.rm = TRUE)
    pval <- stats::wilcox.test(g0, g1)$p.value
    tibble(
      Variable = var_name,
      Overall_N_pct_or_Median_IQR = sprintf("%.1f [%.1f - %.1f]", mt, qt_lo, qt_hi),
      CECI_Low_N_pct_or_Median_IQR = sprintf("%.1f [%.1f - %.1f]", m0, q0_lo, q0_hi),
      CECI_High_N_pct_or_Median_IQR = sprintf("%.1f [%.1f - %.1f]", m1, q1_lo, q1_hi),
      P_Value = pval
    )
  }
}

t1_vars <- list(
  c("age", FALSE),
  c("sex", TRUE), # Handled as male=1
  c("vasopressor_24h", TRUE),
  c("imv_24h", TRUE),
  c("rrt_24h", TRUE),
  c("shock_flag", TRUE),
  c("lactate_24h", FALSE),
  c("albumin_24h", FALSE),
  c("platelet_24h", FALSE),
  c("inr_24h", FALSE),
  c("lar_24h", FALSE),
  c("creatinine_24h", FALSE),
  c("map_mean_24h", FALSE),
  c("charlson_approx", FALSE),
  c("cirrhosis_portal_htn", TRUE),
  c("ckd", TRUE),
  c("chf", TRUE),
  c("diabetes", TRUE),
  c("malignancy", TRUE),
  c("source_pulmonary", TRUE),
  c("source_abdominal", TRUE),
  c("source_urinary", TRUE),
  c("icu_los_days", FALSE),
  c("mort_28d", TRUE),
  c("mort_90d", TRUE)
)

df_primary <- df_primary %>% dplyr::mutate(sex_male = ifelse(sex == "Male", 1L, 0L))
t1_vars[[2]][1] <- "sex_male"

t1_list <- lapply(t1_vars, function(v) make_table1_row(df_primary, v[1], as.logical(v[2])))
table1_df <- dplyr::bind_rows(t1_list)
readr::write_csv(table1_df, file.path(TAB_DIR, "Table_1_MIMIC_Baseline_Characteristics.csv"))
cat("Table 1 exported successfully.\n")

# ==============================================================================
# Step 3: Multivariable Cox Proportional Hazards Models
# ==============================================================================
cat("\n=== Step 3: Fitting Multivariable Cox Proportional Hazards Models ===\n")

# Model 0: Unadjusted
fit_m0 <- survival::coxph(Surv(time_28d, mort_28d) ~ ceci_z, data = df_primary)
fit_m0_bin <- survival::coxph(Surv(time_28d, mort_28d) ~ ceci_high, data = df_primary)
fit_m0_q <- survival::coxph(Surv(time_28d, mort_28d) ~ ceci_quartile, data = df_primary)

# Model 1: Demographics (Age, Sex, Race)
fit_m1 <- survival::coxph(
  Surv(time_28d, mort_28d) ~ ceci_z + age + sex + race_clean,
  data = df_primary
)

# Model 2: Demographics + Acuity (Vasopressor, IMV, RRT, Care Unit)
fit_m2 <- survival::coxph(
  Surv(time_28d, mort_28d) ~ ceci_z + age + sex + race_clean + vasopressor_24h + imv_24h + rrt_24h + careunit_group,
  data = df_primary
)

# Model 3: Fully Adjusted (Acuity + Comorbidities + Infection Source + Other Organ Labs)
fit_m3 <- survival::coxph(
  Surv(time_28d, mort_28d) ~ ceci_z + age + sex + race_clean + vasopressor_24h + imv_24h + rrt_24h + careunit_group +
    charlson_approx + cirrhosis_portal_htn + ckd + chf + diabetes + malignancy +
    source_pulmonary + source_abdominal + source_urinary +
    creatinine_24h + map_mean_imp + bilirubin_imp + bilirubin_na + wbc_imp,
  data = df_primary
)

fit_m3_bin <- survival::coxph(
  Surv(time_28d, mort_28d) ~ ceci_high + age + sex + race_clean + vasopressor_24h + imv_24h + rrt_24h + careunit_group +
    charlson_approx + cirrhosis_portal_htn + ckd + chf + diabetes + malignancy +
    source_pulmonary + source_abdominal + source_urinary +
    creatinine_24h + map_mean_imp + bilirubin_imp + bilirubin_na + wbc_imp,
  data = df_primary
)

extract_cox_res <- function(fit, term_name, model_label) {
  s <- summary(fit)
  cf <- s$coefficients[term_name, ]
  ci <- s$conf.int[term_name, ]
  tibble(
    Model = model_label,
    Term = term_name,
    N = fit$n,
    Events = fit$nevent,
    HR = unname(ci[1]),
    HR_Lower = unname(ci[3]),
    HR_Upper = unname(ci[4]),
    P_Value = unname(cf[5]),
    Concordance = unname(s$concordance[1]),
    Concordance_SE = unname(s$concordance[2])
  )
}

table2_rows <- list(
  extract_cox_res(fit_m0, "ceci_z", "Model 0 (Unadjusted, continuous per 1-SD)"),
  extract_cox_res(fit_m0_bin, "ceci_high", "Model 0 (Unadjusted, High vs Low CECI)"),
  extract_cox_res(fit_m1, "ceci_z", "Model 1 (+ Demographics)"),
  extract_cox_res(fit_m2, "ceci_z", "Model 2 (+ ICU Care Unit & Life Support)"),
  extract_cox_res(fit_m3, "ceci_z", "Model 3 (Fully Adjusted Multivariable)"),
  extract_cox_res(fit_m3_bin, "ceci_high", "Model 3 (Fully Adjusted, High vs Low CECI)")
)
table2_df <- dplyr::bind_rows(table2_rows)
readr::write_csv(table2_df, file.path(TAB_DIR, "Table_2_MIMIC_Multivariable_Cox_Models.csv"))
cat("Table 2 Multivariable Cox summary exported:\n")
print(table2_df %>% dplyr::select(Model, HR, HR_Lower, HR_Upper, P_Value, Concordance) %>% as.data.frame())

# ==============================================================================
# Step 4: Observational Association via Propensity Score IPTW
# ==============================================================================
cat("\n=== Step 4: Observational Association via Propensity Score IPTW ===\n")

# Fit Propensity Score for Elevated CECI (ceci_high = 1)
ps_formula <- ceci_high ~ age + sex + race_clean + careunit_group +
  vasopressor_24h + imv_24h + rrt_24h + shock_flag +
  charlson_approx + cirrhosis_portal_htn + ckd + chf + diabetes + malignancy +
  source_pulmonary + source_abdominal + source_urinary +
  creatinine_24h + map_mean_imp + bilirubin_imp + bilirubin_na

fit_ps <- stats::glm(ps_formula, data = df_primary, family = stats::binomial)
df_primary$ps <- stats::predict(fit_ps, type = "response")

# Stabilized Inverse Probability of Treatment Weights (SW)
p_exp <- mean(df_primary$ceci_high)
df_primary$sw_raw <- ifelse(
  df_primary$ceci_high == 1,
  p_exp / df_primary$ps,
  (1 - p_exp) / (1 - df_primary$ps)
)

# Truncate extreme weights at 1st and 99th percentiles
q01 <- stats::quantile(df_primary$sw_raw, 0.01)
q99 <- stats::quantile(df_primary$sw_raw, 0.99)
df_primary$sw <- pmin(pmax(df_primary$sw_raw, q01), q99)

cat(sprintf("Propensity Score range: %.4f - %.4f (mean = %.4f)\n", min(df_primary$ps), max(df_primary$ps), mean(df_primary$ps)))
cat(sprintf("Stabilized Weights: raw mean = %.4f, truncated [1%%-99%%] mean = %.4f (sd = %.4f)\n",
            mean(df_primary$sw_raw), mean(df_primary$sw), stats::sd(df_primary$sw)))

# Calculate Standardized Mean Differences (SMD) before and after IPTW
cov_names <- c("age", "vasopressor_24h", "imv_24h", "rrt_24h", "charlson_approx",
               "cirrhosis_portal_htn", "ckd", "chf", "diabetes", "malignancy",
               "source_pulmonary", "source_abdominal", "source_urinary",
               "creatinine_24h", "map_mean_imp")

calc_smd <- function(df, var, wt_col = NULL) {
  if (is.null(wt_col)) {
    x0 <- df %>% dplyr::filter(ceci_high == 0) %>% dplyr::pull(!!sym(var))
    x1 <- df %>% dplyr::filter(ceci_high == 1) %>% dplyr::pull(!!sym(var))
    m0 <- mean(x0, na.rm = TRUE); v0 <- stats::var(x0, na.rm = TRUE)
    m1 <- mean(x1, na.rm = TRUE); v1 <- stats::var(x1, na.rm = TRUE)
    s_pool <- sqrt((v0 + v1) / 2)
    abs(m1 - m0) / s_pool
  } else {
    w <- df[[wt_col]]
    df0 <- df %>% dplyr::filter(ceci_high == 0); w0 <- df0[[wt_col]]
    df1 <- df %>% dplyr::filter(ceci_high == 1); w1 <- df1[[wt_col]]
    m0 <- sum(w0 * df0[[var]], na.rm = TRUE) / sum(w0)
    m1 <- sum(w1 * df1[[var]], na.rm = TRUE) / sum(w1)
    v0 <- sum(w0 * (df0[[var]] - m0)^2, na.rm = TRUE) / sum(w0)
    v1 <- sum(w1 * (df1[[var]] - m1)^2, na.rm = TRUE) / sum(w1)
    s_pool <- sqrt((v0 + v1) / 2)
    abs(m1 - m0) / s_pool
  }
}

smd_df <- tibble(
  Covariate = cov_names,
  Unweighted_SMD = vapply(cov_names, function(v) calc_smd(df_primary, v), numeric(1)),
  Weighted_IPTW_SMD = vapply(cov_names, function(v) calc_smd(df_primary, v, "sw"), numeric(1))
)
readr::write_csv(smd_df, file.path(TAB_DIR, "Table_S2_MIMIC_IPTW_Covariate_Balance.csv"))
cat("Covariate balance audited: Max unweighted SMD = ", max(smd_df$Unweighted_SMD), " -> Max IPTW SMD = ", max(smd_df$Weighted_IPTW_SMD), "\n")

# Fit IPTW-Weighted Cox Model with Robust Sandwich SE
fit_iptw <- survival::coxph(
  Surv(time_28d, mort_28d) ~ ceci_high,
  data = df_primary,
  weights = df_primary$sw,
  robust = TRUE
)
s_iptw <- summary(fit_iptw)
iptw_row <- tibble(
  Method = "Propensity Score Stabilized IPTW (Truncated 1-99%)",
  Exposure = "High CECI (>= 0) vs Low CECI (< 0)",
  N = fit_iptw$n,
  Weighted_N = sum(df_primary$sw),
  Events = fit_iptw$nevent,
  HR = unname(s_iptw$conf.int[1, 1]),
  HR_Lower = unname(s_iptw$conf.int[1, 3]),
  HR_Upper = unname(s_iptw$conf.int[1, 4]),
  Robust_SE = unname(s_iptw$coefficients[1, 4]),
  P_Value = unname(s_iptw$coefficients[1, 6])
)
readr::write_csv(iptw_row, file.path(TAB_DIR, "Table_3_MIMIC_IPTW_Association.csv"))
cat("Table 3 IPTW Association exported:\n")
print(iptw_row %>% as.data.frame())

# ==============================================================================
# Step 5: Sensitivity Analyses (Landmark, Subgroups, Extended Cohort, 90d)
# ==============================================================================
cat("\n=== Step 5: Sensitivity & Subgroup Analyses ===\n")

# 1. Confirm survival from hour 24 within the already eligible 24h-survivor cohort
df_landmark <- df_primary %>% dplyr::filter(time_28d > 1.0)
fit_landmark <- survival::coxph(
  Surv(time_post24h_27d, mort_28d) ~ ceci_z + age + sex + race_clean + vasopressor_24h + imv_24h + rrt_24h + careunit_group +
    charlson_approx + cirrhosis_portal_htn + ckd + chf + diabetes + malignancy +
    source_pulmonary + source_abdominal + source_urinary +
    creatinine_24h + map_mean_imp + bilirubin_imp + bilirubin_na + wbc_imp,
  data = df_landmark
)

# 2. 90-Day Mortality
fit_90d <- survival::coxph(
  Surv(time_90d, mort_90d) ~ ceci_z + age + sex + race_clean + vasopressor_24h + imv_24h + rrt_24h + careunit_group +
    charlson_approx + cirrhosis_portal_htn + ckd + chf + diabetes + malignancy +
    source_pulmonary + source_abdominal + source_urinary +
    creatinine_24h + map_mean_imp + bilirubin_imp + bilirubin_na + wbc_imp,
  data = df_primary
)

# 3. Subgroup Interactions
subgroups <- list(
  list(name = "Age >= 65", subset = df_primary$age >= 65),
  list(name = "Age < 65",  subset = df_primary$age < 65),
  list(name = "Vasopressor use", subset = df_primary$vasopressor_24h == 1),
  list(name = "No Vasopressor", subset = df_primary$vasopressor_24h == 0),
  list(name = "Invasive Mechanical Ventilation", subset = df_primary$imv_24h == 1),
  list(name = "No Mechanical Ventilation", subset = df_primary$imv_24h == 0),
  list(name = "With Cirrhosis/Liver Disease", subset = df_primary$cirrhosis_portal_htn == 1),
  list(name = "Without Cirrhosis", subset = df_primary$cirrhosis_portal_htn == 0),
  list(name = "Pulmonary Sepsis", subset = df_primary$source_pulmonary == 1),
  list(name = "Abdominal Sepsis", subset = df_primary$source_abdominal == 1)
)

sub_rows <- list()
for (sg in subgroups) {
  sub_data <- df_primary[sg$subset, ]
  fit_sg <- survival::coxph(Surv(time_28d, mort_28d) ~ ceci_z + age + sex + vasopressor_24h + imv_24h, data = sub_data)
  s_sg <- summary(fit_sg)
  sub_rows[[length(sub_rows) + 1L]] <- tibble(
    Subgroup = sg$name,
    N = nrow(sub_data),
    Events = sum(sub_data$mort_28d),
    HR = unname(s_sg$conf.int["ceci_z", 1]),
    HR_Lower = unname(s_sg$conf.int["ceci_z", 3]),
    HR_Upper = unname(s_sg$conf.int["ceci_z", 4]),
    P_Value = unname(s_sg$coefficients["ceci_z", 5])
  )
}
sub_df <- dplyr::bind_rows(sub_rows)
readr::write_csv(sub_df, file.path(TAB_DIR, "Table_S3_MIMIC_Subgroup_Analyses.csv"))

# Sensitivity Summary Table
sens_df <- tibble(
  Sensitivity_Analysis = c(
    sprintf("Primary 28-day Complete Panel (N=%d, %d deaths)", fit_m3$n, fit_m3$nevent),
    sprintf("Post-24h time origin (same already eligible population; N=%d, %d deaths)", fit_landmark$n, fit_landmark$nevent),
    sprintf("90-Day Extended Follow-up (N=%d, %d deaths)", fit_90d$n, fit_90d$nevent),
    "Propensity Score Weighted High-vs-Low CECI Association"
  ),
  HR = c(
    unname(summary(fit_m3)$conf.int["ceci_z", 1]),
    unname(summary(fit_landmark)$conf.int["ceci_z", 1]),
    unname(summary(fit_90d)$conf.int["ceci_z", 1]),
    iptw_row$HR
  ),
  HR_Lower = c(
    unname(summary(fit_m3)$conf.int["ceci_z", 3]),
    unname(summary(fit_landmark)$conf.int["ceci_z", 3]),
    unname(summary(fit_90d)$conf.int["ceci_z", 3]),
    iptw_row$HR_Lower
  ),
  HR_Upper = c(
    unname(summary(fit_m3)$conf.int["ceci_z", 4]),
    unname(summary(fit_landmark)$conf.int["ceci_z", 4]),
    unname(summary(fit_90d)$conf.int["ceci_z", 4]),
    iptw_row$HR_Upper
  ),
  P_Value = c(
    unname(summary(fit_m3)$coefficients["ceci_z", 5]),
    unname(summary(fit_landmark)$coefficients["ceci_z", 5]),
    unname(summary(fit_90d)$coefficients["ceci_z", 5]),
    iptw_row$P_Value
  )
)
readr::write_csv(sens_df, file.path(TAB_DIR, "Table_S4_MIMIC_Sensitivity_Summary.csv"))
cat("Table S4 Sensitivity Summary exported.\n")

# ==============================================================================
# Step 6: Generating Publication-Grade Figures (Fig 3, Fig 4, Fig 5)
# ==============================================================================
cat("\n=== Step 6: Rendering Publication Figures (Fig 3 - Fig 5) ===\n")

theme_nature <- function() {
  theme_classic(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold", size = 12, hjust = 0),
      plot.subtitle = element_text(size = 10, color = "grey30"),
      axis.title = element_text(face = "bold", size = 11),
      axis.text = element_text(size = 9, color = "black"),
      legend.title = element_text(face = "bold", size = 10),
      legend.text = element_text(size = 9),
      panel.grid.major = element_line(color = "grey92", linewidth = 0.3),
      panel.grid.minor = element_blank()
    )
}

# --- Figure 3: MIMIC-IV Clinical Landscape & IPTW Balance ---
# 3A: Propensity score overlap
p3a <- ggplot(df_primary, aes(x = ps, fill = factor(ceci_high, labels = c("Low CECI (< 0)", "High CECI (>= 0)")))) +
  geom_density(alpha = 0.5, color = "black", linewidth = 0.4) +
  scale_fill_manual(values = c("#2166ac", "#b2182b"), name = "Exposure") +
  labs(title = "A. Propensity Score Overlap", x = "Estimated Propensity Score", y = "Density") +
  theme_nature() +
  theme(legend.position = c(0.75, 0.8))

# 3B: Love Plot of Covariate Balance before and after IPTW
love_data <- smd_df %>%
  tidyr::pivot_longer(cols = c("Unweighted_SMD", "Weighted_IPTW_SMD"), names_to = "Status", values_to = "SMD") %>%
  dplyr::mutate(Status = ifelse(Status == "Unweighted_SMD", "Unadjusted", "IPTW Weighted"))

p3b <- ggplot(love_data, aes(x = SMD, y = reorder(Covariate, SMD), color = Status, shape = Status)) +
  geom_vline(xintercept = 0.10, linetype = "dashed", color = "firebrick", linewidth = 0.5) +
  geom_point(size = 2.5) +
  scale_color_manual(values = c("Unadjusted" = "#b2182b", "IPTW Weighted" = "#2166ac")) +
  scale_shape_manual(values = c(16, 17)) +
  labs(title = "B. Love Plot (Covariate Balance)", x = "Absolute Standardized Mean Difference", y = "") +
  theme_nature() +
  theme(legend.position = c(0.75, 0.2))

# 3C: CECI Z-score distribution by 28-day survival
p3c <- ggplot(df_primary, aes(x = factor(mort_28d, labels = c("28-Day Survivor", "28-Day Non-Survivor")), y = ceci_z, fill = factor(mort_28d))) +
  geom_boxplot(width = 0.4, alpha = 0.7, outlier.shape = 16, outlier.size = 1) +
  scale_fill_manual(values = c("#4393c3", "#d6604d"), guide = "none") +
  labs(title = "C. Clinical Coagulopathy Index", x = "", y = "CECI Z-Score (Day-1)") +
  theme_nature()

# 3D: Correlation Matrix of Panel Biomarkers
biomarkers <- df_primary %>% dplyr::select(lactate_24h, albumin_24h, platelet_24h, inr_24h, ceci_z)
cor_mat <- cor(biomarkers, use = "complete.obs")
cor_df <- as.data.frame(as.table(cor_mat))
colnames(cor_df) <- c("Var1", "Var2", "Corr")

p3d <- ggplot(cor_df, aes(x = Var1, y = Var2, fill = Corr)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%.2f", Corr)), size = 3) +
  scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b", midpoint = 0, limit = c(-1, 1)) +
  labs(title = "D. Biomarker Correlation", x = "", y = "") +
  theme_nature() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

fig3 <- (p3a | p3b) / (p3c | p3d)
ggsave(file.path(FIG_DIR, "Fig3_MIMIC_Clinical_Landscape_and_IPTW.pdf"), fig3, width = 12, height = 10)
ggsave(file.path(FIG_DIR, "Fig3_MIMIC_Clinical_Landscape_and_IPTW.png"), fig3, width = 12, height = 10, dpi = 300)
cat("Fig3 saved successfully.\n")

# --- Figure 4: Kaplan-Meier Survival Curves & Multivariable Cox Forest Plot ---
km_fit_ceci <- survival::survfit(Surv(time_28d, mort_28d) ~ ceci_high, data = df_primary)
km_test <- survival::survdiff(Surv(time_28d, mort_28d) ~ ceci_high, data = df_primary)
km_p <- stats::pchisq(km_test$chisq, df = 1, lower.tail = FALSE)
km_df <- tibble(
  time = km_fit_ceci$time,
  surv = km_fit_ceci$surv,
  lower = km_fit_ceci$lower,
  upper = km_fit_ceci$upper,
  strata = rep(c("Low CECI (< 0)", "High CECI (>= 0)"), km_fit_ceci$strata)
)

p4a <- ggplot(km_df, aes(x = time, y = surv, color = strata, fill = strata)) +
  geom_step(linewidth = 0.8) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.15, color = NA) +
  scale_color_manual(values = c("Low CECI (< 0)" = "#2166ac", "High CECI (>= 0)" = "#b2182b"), name = "CECI Phenotype") +
  scale_fill_manual(values = c("Low CECI (< 0)" = "#2166ac", "High CECI (>= 0)" = "#b2182b"), name = "CECI Phenotype") +
  scale_y_continuous(limits = c(0.4, 1.0), labels = scales::percent) +
  labs(title = "A. Kaplan-Meier 28-Day Survival in MIMIC-IV (N = 3,254)",
       subtitle = sprintf("24h-survivor cohort | Log-rank P = %.2e", km_p),
       x = "Days from ICU Admission", y = "Survival Probability") +
  theme_nature() +
  theme(legend.position = c(0.25, 0.25))

# Forest plot of multivariable models
forest_data <- table2_df %>%
  dplyr::filter(Term %in% c("ceci_z", "ceci_high")) %>%
  dplyr::mutate(
    Model_Clean = c(
      "Unadjusted per 1-SD",
      "Unadjusted High vs Low",
      "Model 1 (+ Demographics)",
      "Model 2 (+ ICU Life Support)",
      "Model 3 (Fully Adjusted)",
      "Model 3 High vs Low"
    ),
    order_id = dplyr::row_number()
  )

p4b <- ggplot(forest_data, aes(x = HR, y = reorder(Model_Clean, -order_id), xmin = HR_Lower, xmax = HR_Upper)) +
  geom_vline(xintercept = 1.0, linetype = "dashed", color = "grey50") +
  geom_errorbar(aes(xmin = HR_Lower, xmax = HR_Upper), width = 0.25, linewidth = 0.6, color = "navy") +
  geom_point(size = 3.2, color = "navy") +
  geom_text(aes(label = sprintf("%.2f (%.2f-%.2f)", HR, HR_Lower, HR_Upper)), vjust = -1.2, size = 3) +
  scale_x_log10(breaks = c(1.0, 1.5, 2.0, 2.5, 3.0, 4.0)) +
  labs(title = "B. Multivariable Cox Hazard Ratios for 28-Day Mortality",
       subtitle = "Adjustment sets include demographics, organ support, laboratories, and comorbidities",
       x = "Hazard Ratio (95% CI)", y = "") +
  theme_nature()

fig4 <- p4a / p4b
ggsave(file.path(FIG_DIR, "Fig4_MIMIC_Survival_Curves_and_Cox_Forest.pdf"), fig4, width = 11, height = 10)
ggsave(file.path(FIG_DIR, "Fig4_MIMIC_Survival_Curves_and_Cox_Forest.png"), fig4, width = 11, height = 10, dpi = 300)
cat("Fig4 saved successfully.\n")

# --- Figure 5: Restricted Cubic Splines & Subgroup Interaction Forest Plot ---
# Fit restricted cubic splines for non-linear hazard
fit_rcs <- survival::coxph(
  Surv(time_28d, mort_28d) ~ splines::ns(ceci_z, df = 3) + age + sex + vasopressor_24h + imv_24h,
  data = df_primary
)

ceci_grid <- seq(stats::quantile(df_primary$ceci_z, 0.01), stats::quantile(df_primary$ceci_z, 0.99), length.out = 100)
reference <- data.frame(ceci_z = 0, age = 65, sex = "Male", vasopressor_24h = 0, imv_24h = 0)
prediction_data <- data.frame(ceci_z = ceci_grid, age = 65, sex = "Male", vasopressor_24h = 0, imv_24h = 0)
design_delta <- sweep(model.matrix(fit_rcs, data = prediction_data), 2,
                      as.numeric(model.matrix(fit_rcs, data = reference)), "-")
log_hr <- as.numeric(design_delta %*% stats::coef(fit_rcs))
contrast_se <- sqrt(rowSums((design_delta %*% stats::vcov(fit_rcs)) * design_delta))
fit_linear <- update(fit_rcs, . ~ . - splines::ns(ceci_z, df = 3) + ceci_z)
fit_no_ceci <- update(fit_rcs, . ~ . - splines::ns(ceci_z, df = 3))
p_overall <- stats::pchisq(2 * (as.numeric(logLik(fit_rcs)) - as.numeric(logLik(fit_no_ceci))), 3, lower.tail = FALSE)
p_nonlinear <- stats::pchisq(2 * (as.numeric(logLik(fit_rcs)) - as.numeric(logLik(fit_linear))), 2, lower.tail = FALSE)
readr::write_csv(tibble(Test = c("CECI overall natural cubic spline (3 df)", "Nonlinearity versus linear CECI (2 df)", "High-versus-low CECI log-rank"),
                        P_Value = c(p_overall, p_nonlinear, km_p)), file.path(TAB_DIR, "Table_S10_CECI_Spline_and_Logrank_Tests.csv"))

spline_df <- tibble(
  ceci_z = ceci_grid,
  log_hr = log_hr,
  hr = exp(log_hr),
  hr_lower = exp(log_hr - 1.96 * contrast_se),
  hr_upper = exp(log_hr + 1.96 * contrast_se)
)

p5a <- ggplot(spline_df, aes(x = ceci_z, y = hr)) +
  geom_hline(yintercept = 1.0, linetype = "dashed", color = "grey50") +
  geom_ribbon(aes(ymin = hr_lower, ymax = hr_upper), alpha = 0.2, fill = "firebrick") +
  geom_line(color = "firebrick", linewidth = 0.9) +
  scale_y_log10() +
  labs(title = "A. CECI and Mortality (Natural Cubic Spline, 3 df)",
       subtitle = sprintf("Reference: CECI = 0 | Overall P = %.2e | Nonlinearity P = %.3g", p_overall, p_nonlinear),
       x = "Clinical Coagulopathy Index (CECI Z-Score)", y = "Hazard Ratio for 28-Day Mortality") +
  theme_nature()

# Subgroup Forest Plot
sub_plot_df <- sub_df %>% dplyr::mutate(order_id = dplyr::row_number())
p5b <- ggplot(sub_plot_df, aes(x = HR, y = reorder(Subgroup, -order_id), xmin = HR_Lower, xmax = HR_Upper)) +
  geom_vline(xintercept = 1.0, linetype = "dashed", color = "grey50") +
  geom_errorbar(aes(xmin = HR_Lower, xmax = HR_Upper), width = 0.25, linewidth = 0.5, color = "#2166ac") +
  geom_point(size = 2.8, color = "#2166ac") +
  geom_text(aes(label = sprintf("%.2f (%.2f-%.2f)", HR, HR_Lower, HR_Upper)), vjust = -1.2, size = 2.7) +
  scale_x_log10(breaks = c(1.0, 1.5, 2.0, 2.5, 3.0)) +
  labs(title = "B. Stratum-specific CECI associations",
       subtitle = "Stratum-specific estimates; no formal interaction test in this panel",
       x = "Adjusted Hazard Ratio per 1-SD CECI", y = "") +
  theme_nature()

fig5 <- p5a / p5b
ggsave(file.path(FIG_DIR, "Fig5_MIMIC_Splines_and_Subgroups.pdf"), fig5, width = 11, height = 10)
ggsave(file.path(FIG_DIR, "Fig5_MIMIC_Splines_and_Subgroups.png"), fig5, width = 11, height = 10, dpi = 300)
cat("Fig5 saved successfully.\n")

cat("\n=== 04_mimic_cox_iptw_and_figures.R Completed 100% Successfully ===\n")
