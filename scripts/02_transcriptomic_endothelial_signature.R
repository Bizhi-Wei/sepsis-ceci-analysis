#!/usr/bin/env Rscript
# ==============================================================================
# 02_transcriptomic_endothelial_signature.R
# Genome-wide differential expression (limma), pathway enrichment (fgsea on
# Hallmark + Reactome), derivation of the locked Endothelial-Immunothrombotic
# Mortality Risk Score (EIMRS) via LASSO-Cox in GSE65682 (n = 479 ICU sepsis,
# 114 deaths) with 1,000x bootstrap optimism correction, multivariable clinical/
# MARS-endotype independence testing, in-cohort thrombocytopenia validation,
# and locked-weight external projection in Day-1 GSE95233 septic shock (n = 51).
# Generates publication-grade Fig1 & Fig2 (.pdf + 300 DPI .png) and CSV tables.
# ==============================================================================

suppressPackageStartupMessages({
  library(limma)
  library(fgsea)
  library(survival)
  library(glmnet)
  library(pROC)
  library(boot)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(tibble)
  library(ggplot2)
  library(ggrepel)
  library(patchwork)
})

set.seed(20261026)

if (!nzchar(Sys.getenv("CECI_WORK_DIR"))) stop("Set CECI_WORK_DIR to a private working folder.")
proj_dir <- normalizePath(Sys.getenv("CECI_WORK_DIR"), mustWork = TRUE)
proc_dir <- file.path(proj_dir, "01_data/processed")
tab_dir  <- file.path(proj_dir, "03_results/tables")
fig_dir  <- file.path(proj_dir, "03_results/figures")
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

cat("[1/7] Loading annotated GEO cohorts and GMT pathway databases...\n")
geo <- readRDS(file.path(proc_dir, "geo_annotated_cohorts.rds"))
shared_genes <- geo$shared_genes
pheno65 <- geo$pheno65
pheno95 <- geo$pheno95
expr65  <- geo$expr65_gene[shared_genes, , drop = FALSE]
expr95  <- geo$expr95_gene[shared_genes, , drop = FALSE]

hallmark_gmt <- fgsea::gmtPathways(file.path(proj_dir, "01_data/pathways/h.all.symbols.gmt"))
reactome_gmt <- suppressWarnings(fgsea::gmtPathways(file.path(proj_dir, "01_data/pathways/reactome2022.gmt")))

# ------------------------------------------------------------------------------
# 2. Genome-wide Differential Expression in GSE65682 (limma empirical Bayes)
# ------------------------------------------------------------------------------
cat("[2/7] Running limma differential expression in GSE65682...\n")

# 2A: ICU Sepsis (n=479) vs Healthy Controls (n=42)
idx_sh <- which(pheno65$cohort_group %in% c("Healthy_Control", "Sepsis_28d_Prognostic"))
grp_sh <- factor(
  ifelse(pheno65$cohort_group[idx_sh] == "Sepsis_28d_Prognostic", "Sepsis", "Healthy"),
  levels = c("Healthy", "Sepsis")
)
fit_sh <- limma::eBayes(limma::lmFit(expr65[, idx_sh, drop = FALSE], model.matrix(~ grp_sh)))
deg_sepsis_vs_healthy <- limma::topTable(fit_sh, coef = 2, number = Inf, sort.by = "none") %>%
  tibble::rownames_to_column("gene") %>%
  dplyr::arrange(adj.P.Val, P.Value)

# 2B: ICU Sepsis (n=479) vs Non-Infectious GI Surgical ICU Controls (n=42)
idx_sni <- which(pheno65$cohort_group %in% c("NonInfectious_ICU_Control", "Sepsis_28d_Prognostic"))
grp_sni <- factor(
  ifelse(pheno65$cohort_group[idx_sni] == "Sepsis_28d_Prognostic", "Sepsis", "NonInf_ICU"),
  levels = c("NonInf_ICU", "Sepsis")
)
fit_sni <- limma::eBayes(limma::lmFit(expr65[, idx_sni, drop = FALSE], model.matrix(~ grp_sni)))
deg_sepsis_vs_noninf <- limma::topTable(fit_sni, coef = 2, number = Inf, sort.by = "none") %>%
  tibble::rownames_to_column("gene") %>%
  dplyr::arrange(adj.P.Val, P.Value)

# 2C: 28-Day Non-Survivors (n=114) vs Survivors (n=365) within GSE65682 Sepsis (n=479)
idx_sep <- which(pheno65$cohort_group == "Sepsis_28d_Prognostic")
sepsis_meta65 <- pheno65[idx_sep, ]
sepsis_expr65 <- expr65[, idx_sep, drop = FALSE]

grp_mort <- factor(
  ifelse(sepsis_meta65$mort_28d == 1, "NonSurvivor", "Survivor"),
  levels = c("Survivor", "NonSurvivor")
)
fit_mort <- limma::eBayes(limma::lmFit(sepsis_expr65, model.matrix(~ grp_mort)))
deg_nonsurv_vs_surv <- limma::topTable(fit_mort, coef = 2, number = Inf, sort.by = "none") %>%
  tibble::rownames_to_column("gene") %>%
  dplyr::arrange(P.Value)

readr::write_csv(deg_sepsis_vs_healthy, file.path(tab_dir, "Table_S5_GSE65682_DEG_Sepsis_vs_Healthy.csv"))
readr::write_csv(deg_nonsurv_vs_surv, file.path(tab_dir, "Table_S6_GSE65682_DEG_28dNonSurvivor_vs_Survivor.csv"))

cat(sprintf("  Sepsis vs Healthy DEGs (FDR < 0.05 & |log2FC| > 0.5): %d\n",
            sum(deg_sepsis_vs_healthy$adj.P.Val < 0.05 & abs(deg_sepsis_vs_healthy$logFC) > 0.5)))
cat(sprintf("  28d Non-Survivor vs Survivor DEGs (FDR < 0.05): %d | (FDR < 0.10): %d | (P < 0.01): %d\n",
            sum(deg_nonsurv_vs_surv$adj.P.Val < 0.05),
            sum(deg_nonsurv_vs_surv$adj.P.Val < 0.10),
            sum(deg_nonsurv_vs_surv$P.Value < 0.01)))

# ------------------------------------------------------------------------------
# 3. Genome-wide Pathway Enrichment (fgsea: Hallmark + Reactome)
# ------------------------------------------------------------------------------
cat("[3/7] Running fgsea pathway enrichment on 28-day mortality & sepsis ranks...\n")
ranks_mort <- deg_nonsurv_vs_surv$t
names(ranks_mort) <- deg_nonsurv_vs_surv$gene
ranks_mort <- sort(ranks_mort, decreasing = TRUE)

ranks_sh <- deg_sepsis_vs_healthy$t
names(ranks_sh) <- deg_sepsis_vs_healthy$gene
ranks_sh <- sort(ranks_sh, decreasing = TRUE)

if (Sys.getenv("EIMRS_REUSE_ENRICHMENT", "0") == "1") {
  fg_all <- readr::read_csv(file.path(tab_dir, "Table_S7_FGSEA_Pathway_Enrichment.csv"), show_col_types = FALSE)
  fg_hall_mort <- dplyr::filter(fg_all, database == "Hallmark", comparison == "28d_NonSurvivor_vs_Survivor")
  fg_reac_mort <- dplyr::filter(fg_all, database == "Reactome", comparison == "28d_NonSurvivor_vs_Survivor")
  fg_hall_sh <- dplyr::filter(fg_all, database == "Hallmark", comparison == "Sepsis_vs_Healthy")
  cat("  Reusing retained enrichment output for figure-only revision.\n")
} else {
fg_hall_mort <- fgsea::fgseaMultilevel(pathways = hallmark_gmt, stats = ranks_mort, minSize = 15, maxSize = 500, eps = 0) %>%
  dplyr::as_tibble() %>%
  dplyr::mutate(database = "Hallmark", comparison = "28d_NonSurvivor_vs_Survivor",
                leadingEdge = vapply(leadingEdge, function(x) paste(head(x, 12), collapse = ";"), character(1))) %>%
  dplyr::arrange(padj, pval)

fg_reac_mort <- fgsea::fgseaMultilevel(pathways = reactome_gmt, stats = ranks_mort, minSize = 15, maxSize = 500, eps = 0) %>%
  dplyr::as_tibble() %>%
  dplyr::mutate(database = "Reactome", comparison = "28d_NonSurvivor_vs_Survivor",
                leadingEdge = vapply(leadingEdge, function(x) paste(head(x, 12), collapse = ";"), character(1))) %>%
  dplyr::arrange(padj, pval)

fg_hall_sh <- fgsea::fgseaMultilevel(pathways = hallmark_gmt, stats = ranks_sh, minSize = 15, maxSize = 500, eps = 0) %>%
  dplyr::as_tibble() %>%
  dplyr::mutate(database = "Hallmark", comparison = "Sepsis_vs_Healthy",
                leadingEdge = vapply(leadingEdge, function(x) paste(head(x, 12), collapse = ";"), character(1))) %>%
  dplyr::arrange(padj, pval)

fg_all <- dplyr::bind_rows(fg_hall_mort, fg_reac_mort, fg_hall_sh)
readr::write_csv(fg_all, file.path(tab_dir, "Table_S7_FGSEA_Pathway_Enrichment.csv"))
}

cat("  Top Hallmark pathways in 28d Non-Survivors vs Survivors:\n")
print(fg_hall_mort %>% dplyr::select(pathway, NES, pval, padj) %>% head(10))

# ------------------------------------------------------------------------------
# 4. Curate A Priori Endothelial-Glycocalyx-Immunothrombotic-Hypoxia (EGIH) Universe
#    & Derive Locked LASSO-Cox EIMRS Strictly Inside GSE65682 (n = 479)
# ------------------------------------------------------------------------------
cat("[4/7] Constructing a priori EGIH pathway universe and fitting LASSO-Cox in GSE65682...\n")

egih_hallmark_names <- c(
  "HALLMARK_HYPOXIA",
  "HALLMARK_GLYCOLYSIS",
  "HALLMARK_COAGULATION",
  "HALLMARK_ANGIOGENESIS"
)
egih_reactome_names <- c(
  "Cell Surface Interactions At Vascular Wall R-HSA-202733",
  "VEGFR2 Mediated Vascular Permeability R-HSA-5218920",
  "ECM Proteoglycans R-HSA-3000178",
  "Degradation Of Extracellular Matrix R-HSA-1474228",
  "Formation Of Fibrin Clot (Clotting Cascade) R-HSA-140877",
  "Common Pathway Of Fibrin Clot Formation R-HSA-140875",
  "Extrinsic Pathway Of Fibrin Clot Formation R-HSA-140834",
  "Platelet Degranulation R-HSA-114608",
  "Cellular Response To Hypoxia R-HSA-1234174",
  "Glycolysis R-HSA-70171"
)
curated_vascular_panel <- c(
  "THBD", "TFPI", "PLAT", "SERPINE1", "ANGPT1", "VWF", "PECAM1", "ADM",
  "MMP8", "MMP9", "TIMP1", "S100A8", "S100A12", "IL1R2", "OLFM4", "LCN2",
  "ARG1", "PCSK9", "ITGAM", "ITGAL", "SELL", "CD177", "ELANE", "MPO", "CTSG", "PRTN3",
  "ADAM17", "F5", "F13A1", "HP", "LDHA", "HIF1A", "SLC2A1", "VEGFA", "ICAM1", "ICAM2", "ICAM3", "SDC4", "TIE1", "THBS1"
)

egih_universe <- unique(c(
  unlist(hallmark_gmt[intersect(egih_hallmark_names, names(hallmark_gmt))], use.names = FALSE),
  unlist(reactome_gmt[intersect(egih_reactome_names, names(reactome_gmt))], use.names = FALSE),
  curated_vascular_panel
))
egih_shared <- sort(intersect(egih_universe, shared_genes))
cat(sprintf("  A priori EGIH pathway gene universe in shared blood transcriptome: %d genes.\n", length(egih_shared)))

# Standardize GSE65682 sepsis expression (gene-wise Z-score within GSE65682 n=479 ONLY)
X65_all <- scale(t(sepsis_expr65[egih_shared, , drop = FALSE]))
surv_obj65 <- survival::Surv(time = sepsis_meta65$time_28d, event = sepsis_meta65$mort_28d)

# Univariate Cox & Sepsis-dysregulation screen strictly within GSE65682
sh_lookup <- deg_sepsis_vs_healthy %>% dplyr::select(gene, logFC_SH = logFC, FDR_SH = adj.P.Val)
ns_lookup <- deg_nonsurv_vs_surv %>% dplyr::select(gene, logFC_NS = logFC, P_NS = P.Value, FDR_NS = adj.P.Val)

uni_cox_list <- lapply(egih_shared, function(g) {
  fit <- survival::coxph(surv_obj65 ~ X65_all[, g])
  s <- summary(fit)
  tibble::tibble(
    gene   = g,
    beta   = s$coefficients[1, "coef"],
    HR     = s$coefficients[1, "exp(coef)"],
    HR_lo  = s$conf.int[1, "lower .95"],
    HR_hi  = s$conf.int[1, "upper .95"],
    P_cox  = s$coefficients[1, "Pr(>|z|)"]
  )
})
uni_cox_df <- dplyr::bind_rows(uni_cox_list) %>%
  dplyr::mutate(FDR_cox = p.adjust(P_cox, method = "BH")) %>%
  dplyr::left_join(sh_lookup, by = "gene") %>%
  dplyr::left_join(ns_lookup, by = "gene") %>%
  dplyr::arrange(P_cox)

readr::write_csv(uni_cox_df, file.path(tab_dir, "Table_S8_EGIH_Univariate_Cox_GSE65682.csv"))

# Candidate genes for LASSO-Cox: sepsis-dysregulated (FDR_SH < 0.05 & |logFC_SH| > 0.30)
# AND 28-day mortality-associated (FDR_cox < 0.05 in GSE65682)
cand_df <- uni_cox_df %>%
  dplyr::filter(FDR_SH < 0.05, abs(logFC_SH) > 0.30, FDR_cox < 0.05)
cat(sprintf("  Candidate EGIH genes passing Sepsis-dysregulated & FDR_cox < 0.05 screen in GSE65682: %d genes.\n",
            nrow(cand_df)))

if (nrow(cand_df) < 15) {
  cand_df <- uni_cox_df %>% dplyr::filter(FDR_SH < 0.05, P_cox < 0.005)
  cat(sprintf("  Fallback P_cox < 0.005 candidate count: %d genes.\n", nrow(cand_df)))
}

X65_cand <- X65_all[, cand_df$gene, drop = FALSE]

# Stratified 10-fold CV LASSO-Cox in GSE65682
set.seed(20261026)
foldid <- sample(rep(1:10, length.out = nrow(X65_cand)))
# Balance events across folds
idx_ev0 <- which(sepsis_meta65$mort_28d == 0)
idx_ev1 <- which(sepsis_meta65$mort_28d == 1)
foldid[idx_ev0] <- sample(rep(1:10, length.out = length(idx_ev0)))
foldid[idx_ev1] <- sample(rep(1:10, length.out = length(idx_ev1)))

cv_lasso <- glmnet::cv.glmnet(
  x = X65_cand,
  y = surv_obj65,
  family = "cox",
  alpha = 1,
  foldid = foldid,
  type.measure = "C"
)

# Check non-zero counts at lambda.1se and lambda.min; ensure parsimonious EPV >= 8 (8-16 genes)
coef_1se <- as.matrix(coef(cv_lasso, s = "lambda.1se"))
coef_min <- as.matrix(coef(cv_lasso, s = "lambda.min"))
nz_1se <- sum(coef_1se[, 1] != 0)
nz_min <- sum(coef_min[, 1] != 0)
cat(sprintf("  LASSO-Cox non-zero genes: lambda.1se = %d, lambda.min = %d\n", nz_1se, nz_min))

if (nz_1se >= 6 && nz_1se <= 18) {
  chosen_lambda <- cv_lasso$lambda.1se
  lambda_label <- "lambda.1se"
} else if (nz_min >= 6 && nz_min <= 18) {
  chosen_lambda <- cv_lasso$lambda.min
  lambda_label <- "lambda.min"
} else {
  # Pick largest lambda along path giving between 10 and 15 genes (strict EPV ~ 10 rule for 114 events)
  nz_path <- colSums(as.matrix(cv_lasso$glmnet.fit$beta) != 0)
  ok_idx <- which(nz_path >= 10 & nz_path <= 15)
  if (length(ok_idx) > 0) {
    # Choose one maximizing CV C-index among 10-15 gene models
    best_i <- ok_idx[which.max(cv_lasso$cvm[ok_idx])]
    chosen_lambda <- cv_lasso$lambda[best_i]
    lambda_label <- sprintf("EPV-optimal CV lambda (%d genes)", nz_path[best_i])
  } else {
    chosen_lambda <- cv_lasso$lambda.min
    lambda_label <- "lambda.min"
  }
}

beta_mat <- as.matrix(coef(cv_lasso, s = chosen_lambda))
sel_genes <- rownames(beta_mat)[beta_mat[, 1] != 0]
sel_betas <- as.numeric(beta_mat[sel_genes, 1])
names(sel_betas) <- sel_genes

sig_table <- tibble::tibble(
  gene = sel_genes,
  lasso_beta_per_SD = sel_betas
) %>%
  dplyr::left_join(uni_cox_df, by = "gene") %>%
  dplyr::arrange(dplyr::desc(lasso_beta_per_SD))

cat(sprintf("  Locked EIMRS Signature (%s): %d genes:\n", lambda_label, nrow(sig_table)))
print(sig_table %>% dplyr::select(gene, lasso_beta_per_SD, HR, P_cox, logFC_SH, logFC_NS))
readr::write_csv(sig_table, file.path(tab_dir, "Table_4_Locked_EIMRS_Signature_Genes.csv"))

# Compute EIMRS linear predictor in GSE65682 (raw and Z-standardized)
eimrs65_raw <- as.numeric(X65_all[, sel_genes, drop = FALSE] %*% sel_betas)
eimrs65_z   <- as.numeric(scale(eimrs65_raw))
sepsis_meta65 <- sepsis_meta65 %>%
  dplyr::mutate(
    eimrs_raw = eimrs65_raw,
    eimrs_z   = eimrs65_z,
    eimrs_group = factor(ifelse(eimrs_z >= median(eimrs_z), "High EIMRS", "Low EIMRS"),
                         levels = c("Low EIMRS", "High EIMRS")),
    eimrs_tertile = factor(
      dplyr::ntile(eimrs_z, 3),
      levels = 1:3,
      labels = c("T1 (Low)", "T2 (Intermediate)", "T3 (High)")
    )
  )

# ------------------------------------------------------------------------------
# 5. Rigorous Internal Validation in GSE65682:
#    (a) Apparent & 1,000x Bootstrap Optimism-Corrected C-index & 28d AUC
#    (b) Out-of-fold Nested CV C-index / AUC
#    (c) MARS source subcohorts (n=263 and n=216) Subcohorts
#    (d) Multivariable Cox Independence (beyond Age, Sex, Source, MARS1-4, Thrombocytopenia)
#    (e) In-cohort Biological Link with Clinical Thrombocytopenia Grades
# ------------------------------------------------------------------------------
cat("[5/7] Performing 1,000x bootstrap optimism check & multivariable Cox in GSE65682...\n")

fit_app65 <- survival::coxph(surv_obj65 ~ eimrs_z, data = sepsis_meta65)
c_app65   <- summary(fit_app65)$concordance[1]
c_se65    <- summary(fit_app65)$concordance[2]
roc65     <- pROC::roc(sepsis_meta65$mort_28d, sepsis_meta65$eimrs_z, direction = "<", quiet = TRUE)
auc65_ci  <- as.numeric(pROC::ci.auc(roc65, method = "delong"))

# Full-development validation is generated by 02b_revalidate_eimrs.R.
# Screening, scaling, tuning, and fitting are all repeated in each training set.
validation_path <- file.path(proc_dir, "eimrs_full_pipeline_validation.rds")
if (!file.exists(validation_path)) stop("Run 02b_revalidate_eimrs.R before exporting validation results")
validation <- readRDS(validation_path)
stopifnot(identical(validation$input_md5,
                    unname(tools::md5sum(file.path(proc_dir, "geo_annotated_cohorts.rds")))))
stopifnot(validation$summary$N == nrow(sepsis_meta65),
          validation$summary$Events == sum(sepsis_meta65$mort_28d),
          validation$summary$Full_model_matches_locked_coefficients,
          setequal(names(validation$full_fit$beta), sel_genes),
          max(abs(validation$full_fit$beta[sel_genes] - sel_betas)) < 1e-8,
          abs(validation$summary$Apparent_C - c_app65) < 1e-8)
oof_pred <- validation$oof
oof_c <- validation$summary$OOF_C
oof_auc <- validation$summary$OOF_AUC
B_boot <- validation$summary$Bootstrap_replicates_valid
c_opt_corr <- validation$summary$Corrected_C
roc_oof <- pROC::roc(sepsis_meta65$mort_28d, oof_pred, direction = "<", quiet = TRUE)

cat(sprintf("  GSE65682 Apparent C-index = %.3f (95%% CI %.3f-%.3f) | 1,000x Bootstrap Optimism-Corrected C-index = %.3f\n",
            c_app65, c_app65 - 1.96 * c_se65, c_app65 + 1.96 * c_se65, c_opt_corr))
cat(sprintf("  GSE65682 28d ROC AUC = %.3f (95%% CI %.3f-%.3f) | Nested 10-fold OOF AUC = %.3f (OOF C-index = %.3f)\n",
            auc65_ci[2], auc65_ci[1], auc65_ci[3], oof_auc, oof_c))

# Multivariable Cox models in GSE65682
sepsis_meta65 <- sepsis_meta65 %>%
  dplyr::mutate(
    mars_factor = factor(mars_endotype, levels = c("Mars3", "Mars1", "Mars2", "Mars4")),
    thrombo_clean = ifelse(is.na(thrombocytopenia), "Unrecorded", thrombocytopenia),
    thrombo_factor = factor(thrombo_clean, levels = c("D_normal", "C_low", "B_medium_low", "A_very_low", "Unrecorded"))
  )

cox65_uni  <- survival::coxph(surv_obj65 ~ eimrs_z, data = sepsis_meta65)
cox65_high <- survival::coxph(surv_obj65 ~ eimrs_group, data = sepsis_meta65)
cox65_mv1  <- survival::coxph(surv_obj65 ~ eimrs_z + age + sex + sepsis_source, data = sepsis_meta65)
cox65_mv2  <- survival::coxph(surv_obj65 ~ eimrs_z + age + sex + sepsis_source + mars_factor + thrombo_factor, data = sepsis_meta65)

extract_cox_row <- function(fit, term_name, model_label) {
  s <- summary(fit)
  rn <- grep(term_name, rownames(s$coefficients))[1]
  tibble::tibble(
    Model = model_label,
    Term  = rownames(s$coefficients)[rn],
    N     = fit$n,
    Events = fit$nevent,
    HR    = s$coefficients[rn, "exp(coef)"],
    CI_Lower = s$conf.int[rn, "lower .95"],
    CI_Upper = s$conf.int[rn, "upper .95"],
    P_Value  = s$coefficients[rn, "Pr(>|z|)"],
    C_Index  = s$concordance[1]
  )
}

gse65_cox_summary <- dplyr::bind_rows(
  extract_cox_row(cox65_uni,  "eimrs_z", "GSE65682 Unadjusted (per 1-SD EIMRS)"),
  extract_cox_row(cox65_high, "eimrs_group", "GSE65682 Unadjusted (High vs Low EIMRS)"),
  extract_cox_row(cox65_mv1,  "eimrs_z", "GSE65682 Adjusted: Age + Sex + Sepsis Source"),
  extract_cox_row(cox65_mv2,  "eimrs_z", "GSE65682 Fully Adjusted: + MARS Endotype (Mars1-4) + Thrombocytopenia")
)
print(gse65_cox_summary)

# Check EIMRS across real clinical Thrombocytopenia grades in GSE65682
thrombo_sub <- sepsis_meta65 %>%
  dplyr::filter(thrombocytopenia %in% c("D_normal", "C_low", "B_medium_low", "A_very_low")) %>%
  dplyr::mutate(
    thrombo_ord = factor(thrombocytopenia,
                         levels = c("D_normal", "C_low", "B_medium_low", "A_very_low"),
                         labels = c("Normal", "Mildly Low", "Moderately Low", "Severely Low"))
  )
kw_thrombo <- kruskal.test(eimrs_z ~ thrombo_ord, data = thrombo_sub)
cat(sprintf("  In-cohort link (EIMRS vs Clinical Thrombocytopenia Severity, n=%d): Kruskal-Wallis P = %.3e\n",
            nrow(thrombo_sub), kw_thrombo$p.value))

# ------------------------------------------------------------------------------
# 6. Locked-weight External Projection in GSE95233 Day-1 Septic Shock (n = 51)
#    + Healthy Volunteer Biological Anchor (n = 22)
# ------------------------------------------------------------------------------
cat("[6/7] Running locked-weight external projection in Day-1 GSE95233 (n = 51)...\n")

idx95_d1 <- which(pheno95$cohort_group == "Septic_Shock_Day1")
meta95_d1 <- pheno95[idx95_d1, ]
expr95_d1 <- expr95[sel_genes, idx95_d1, drop = FALSE]

# Standardize within GSE95233 Day-1 baseline cohort (zero outcome label usage!)
X95_d1 <- scale(t(expr95_d1))
eimrs95_raw <- as.numeric(X95_d1 %*% sel_betas)
eimrs95_z   <- as.numeric(scale(eimrs95_raw))

meta95_d1 <- meta95_d1 %>%
  dplyr::mutate(
    eimrs_raw = eimrs95_raw,
    eimrs_z   = eimrs95_z,
    outcome_label = factor(ifelse(mort_event == 1, "Non-Survivor", "Survivor"),
                           levels = c("Survivor", "Non-Survivor"))
  )

roc95 <- pROC::roc(meta95_d1$mort_event, meta95_d1$eimrs_z, direction = "<", quiet = TRUE)
auc95_ci <- as.numeric(pROC::ci.auc(roc95, method = "delong"))
wilcox95_p <- wilcox.test(eimrs_z ~ mort_event, data = meta95_d1)$p.value

glm95_uni <- glm(mort_event ~ eimrs_z, data = meta95_d1, family = binomial())
glm95_adj <- glm(mort_event ~ eimrs_z + age + sex, data = meta95_d1, family = binomial())

s95_uni <- summary(glm95_uni)$coefficients["eimrs_z", ]
s95_adj <- summary(glm95_adj)$coefficients["eimrs_z", ]
or95_uni <- exp(c(s95_uni[1], s95_uni[1] - 1.96 * s95_uni[2], s95_uni[1] + 1.96 * s95_uni[2]))
or95_adj <- exp(c(s95_adj[1], s95_adj[1] - 1.96 * s95_adj[2], s95_adj[1] + 1.96 * s95_adj[2]))

# Project Healthy Controls + Day-1 Survivors + Day-1 Non-Survivors using pooled baseline Z-score for biological gradient
idx95_base <- which(pheno95$cohort_group %in% c("Healthy_Control", "Septic_Shock_Day1"))
meta95_base <- pheno95[idx95_base, ]
X95_base <- scale(t(expr95[sel_genes, idx95_base, drop = FALSE]))
meta95_base$eimrs_base_z <- as.numeric(scale(as.numeric(X95_base %*% sel_betas)))
meta95_base$bio_group <- factor(
  dplyr::case_when(
    meta95_base$cohort_group == "Healthy_Control" ~ "Healthy\n(n=22)",
    meta95_base$mort_event == 0 ~ "Day-1 Shock\nSurvivor (n=34)",
    meta95_base$mort_event == 1 ~ "Day-1 Shock\nNon-Survivor (n=17)"
  ),
  levels = c("Healthy\n(n=22)", "Day-1 Shock\nSurvivor (n=34)", "Day-1 Shock\nNon-Survivor (n=17)")
)
kw95_bio_p <- kruskal.test(eimrs_base_z ~ bio_group, data = meta95_base)$p.value

cat(sprintf("  EXTERNAL VALIDATION (GSE95233 Day-1, n=51, 17 deaths): AUC = %.3f (95%% DeLong CI %.3f-%.3f)\n",
            auc95_ci[2], auc95_ci[1], auc95_ci[3]))
cat(sprintf("  GSE95233 Wilcoxon P = %.4e | Unadjusted OR per 1-SD = %.2f (95%% CI %.2f-%.2f, P = %.4f)\n",
            wilcox95_p, or95_uni[1], or95_uni[2], or95_uni[3], s95_uni[4]))
cat(sprintf("  GSE95233 Age+Sex Adjusted OR per 1-SD = %.2f (95%% CI %.2f-%.2f, P = %.4f)\n",
            or95_adj[1], or95_adj[2], or95_adj[3], s95_adj[4]))

val_summary_df <- tibble::tibble(
  Cohort = c(
    "GSE65682 Discovery (MARS ICU Sepsis, 28d Mortality)",
    sprintf("GSE65682 Full-pipeline Bootstrap Optimism-Corrected (%d replicates)", B_boot),
    "GSE65682 Training-only Screening Nested 10-Fold Out-of-Fold CV",
    "GSE65682 Multivariable Cox (+Age, Sex, Source, MARS1-4, Thrombocytopenia)",
    "GSE95233 Locked-weight External Projection (Day-1, Signature Alone)",
    "GSE95233 Exploratory Age+Sex Augmented Model (Fitted in External Cohort)"
  ),
  Sample_Size_N = c(479, 479, 479, 479, 51, 51),
  Events_n = c(114, 114, 114, 114, 17, 17),
  Discrimination_AUC_or_C = c(auc65_ci[2], c_opt_corr, oof_auc, summary(cox65_mv2)$concordance[1], auc95_ci[2], as.numeric(pROC::auc(pROC::roc(meta95_d1$mort_event, fitted(glm95_adj), quiet = TRUE)))),
  Discrim_CI_Lower = c(auc65_ci[1], NA_real_, NA, summary(cox65_mv2)$concordance[1] - 1.96 * summary(cox65_mv2)$concordance[2], auc95_ci[1], NA),
  Discrim_CI_Upper = c(auc65_ci[3], NA_real_, NA, summary(cox65_mv2)$concordance[1] + 1.96 * summary(cox65_mv2)$concordance[2], auc95_ci[3], NA),
  Effect_Measure = c("HR per 1-SD", "Harrell C-index", "OOF AUC", "Adjusted HR per 1-SD", "OR per 1-SD", "Adjusted OR per 1-SD"),
  Estimate = c(gse65_cox_summary$HR[1], c_opt_corr, oof_auc, gse65_cox_summary$HR[4], or95_uni[1], or95_adj[1]),
  Effect_CI_Lower = c(gse65_cox_summary$CI_Lower[1], NA_real_, NA, gse65_cox_summary$CI_Lower[4], or95_uni[2], or95_adj[2]),
  Effect_CI_Upper = c(gse65_cox_summary$CI_Upper[1], NA_real_, NA, gse65_cox_summary$CI_Upper[4], or95_uni[3], or95_adj[3]),
  P_Value = c(gse65_cox_summary$P_Value[1], NA_real_, NA, gse65_cox_summary$P_Value[4], s95_uni[4], s95_adj[4])
)
readr::write_csv(val_summary_df, file.path(tab_dir, "Table_5_Transcriptomic_EIMRS_Validation_Summary.csv"))
readr::write_csv(gse65_cox_summary, file.path(tab_dir, "Table_S9_GSE65682_Multivariable_Cox_Models.csv"))

# Save scored transcriptomic data for downstream reference
saveRDS(
  list(
    sig_table = sig_table,
    sepsis_meta65 = sepsis_meta65,
    meta95_d1 = meta95_d1,
    meta95_base = meta95_base,
    val_summary_df = val_summary_df,
    gse65_cox_summary = gse65_cox_summary,
    full_pipeline_validation_summary = validation$summary
  ),
  file = file.path(proc_dir, "eimrs_results_bundle.rds")
)

# ------------------------------------------------------------------------------
# 7. Generate Publication-Grade Figures: Fig1 & Fig2 (.pdf + 300 DPI .png)
# ------------------------------------------------------------------------------
cat("[7/7] Rendering publication-grade Fig1 and Fig2...\n")

theme_pub <- function(base_size = 10.5) {
  ggplot2::theme_classic(base_size = base_size) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = base_size + 1, hjust = 0),
      plot.subtitle = ggplot2::element_text(size = base_size - 1, color = "#444444", hjust = 0),
      axis.title = ggplot2::element_text(face = "bold", size = base_size),
      axis.text = ggplot2::element_text(color = "#111111", size = base_size - 1),
      legend.title = ggplot2::element_text(face = "bold", size = base_size - 1),
      legend.text = ggplot2::element_text(size = base_size - 1.5),
      panel.grid.major = ggplot2::element_line(color = "#F0F0F0", linewidth = 0.3)
    )
}

# --- FIG 1: Transcriptomic Endothelial-Immunothrombotic Landscape in Sepsis ---
# Fig 1A: Volcano Plot - Sepsis (n=479) vs Healthy Controls (n=42)
vol_sh <- deg_sepsis_vs_healthy %>%
  dplyr::mutate(
    neglog10FDR = -log10(pmax(adj.P.Val, 1e-300)),
    status = dplyr::case_when(
      adj.P.Val < 0.05 & logFC > 0.5 ~ "Up in Sepsis",
      adj.P.Val < 0.05 & logFC < -0.5 ~ "Down in Sepsis",
      TRUE ~ "Not Significant"
    ),
    label_gene = ifelse(gene %in% unique(c(sel_genes, "MMP8", "OLFM4", "ARG1", "IL1R2", "LCN2", "S100A8", "S100A12", "LDHA", "F5", "TFPI", "PECAM1", "SLC2A1", "ELANE", "MPO")), gene, NA_character_)
  )

p1a <- ggplot2::ggplot(vol_sh, ggplot2::aes(x = logFC, y = neglog10FDR, color = status)) +
  ggplot2::geom_point(alpha = 0.55, size = 1.1) +
  ggplot2::scale_color_manual(values = c("Up in Sepsis" = "#B2182B", "Down in Sepsis" = "#2166AC", "Not Significant" = "#BDBDBD")) +
  ggplot2::geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", color = "#555555", linewidth = 0.4) +
  ggplot2::geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "#555555", linewidth = 0.4) +
  ggrepel::geom_text_repel(
    data = dplyr::filter(vol_sh, !is.na(label_gene)),
    ggplot2::aes(label = label_gene),
    size = 2.8, fontface = "bold.italic", color = "#111111",
    max.overlaps = 22, box.padding = 0.3, segment.size = 0.25
  ) +
  ggplot2::labs(
    title = "Sepsis (n = 479) vs Healthy Controls (n = 42)",
    subtitle = "GSE65682 whole-blood transcriptome (10,817 harmonized genes)",
    x = "log2 Fold Change (Sepsis / Healthy)",
    y = "-log10(Adjusted P-value)",
    color = NULL
  ) +
  theme_pub() +
  ggplot2::theme(legend.position = "top")

# Fig 1B: Prognostic Volcano - 28-Day Non-Survivors (n=114) vs Survivors (n=365)
vol_ns <- deg_nonsurv_vs_surv %>%
  dplyr::mutate(
    neglog10P = -log10(pmax(P.Value, 1e-30)),
    is_egih = gene %in% egih_shared,
    status = dplyr::case_when(
      gene %in% sel_genes & logFC > 0 ~ "Locked EIMRS Risk Gene",
      gene %in% sel_genes & logFC < 0 ~ "Locked EIMRS Protective Gene",
      P.Value < 0.01 & logFC > 0 ~ "Up in Non-Survivors (P < 0.01)",
      P.Value < 0.01 & logFC < 0 ~ "Down in Non-Survivors (P < 0.01)",
      TRUE ~ "Background"
    ),
    label_gene = ifelse(gene %in% unique(c(sel_genes, "SLC2A1", "ELANE", "MPO", "CTSG", "PRTN3", "PECAM1", "F5", "IL1R2", "ARG1", "PCSK9", "VEGFA")), gene, NA_character_)
  )

p1b <- ggplot2::ggplot(vol_ns, ggplot2::aes(x = logFC, y = neglog10P, color = status)) +
  ggplot2::geom_point(alpha = 0.65, size = 1.2) +
  ggplot2::scale_color_manual(values = c(
    "Locked EIMRS Risk Gene" = "#8B0000",
    "Locked EIMRS Protective Gene" = "#08306B",
    "Up in Non-Survivors (P < 0.01)" = "#F4A582",
    "Down in Non-Survivors (P < 0.01)" = "#92C5DE",
    "Background" = "#D9D9D9"
  )) +
  ggplot2::geom_hline(yintercept = -log10(0.01), linetype = "dashed", color = "#555555", linewidth = 0.4) +
  ggplot2::geom_vline(xintercept = 0, linetype = "dotted", color = "#555555", linewidth = 0.4) +
  ggrepel::geom_text_repel(
    data = dplyr::filter(vol_ns, !is.na(label_gene)),
    ggplot2::aes(label = label_gene),
    size = 2.8, fontface = "bold.italic", color = "#111111",
    max.overlaps = 25, box.padding = 0.3, segment.size = 0.25
  ) +
  ggplot2::labs(
    title = "28-Day Non-Survivors (n = 114) vs Survivors (n = 365)",
    subtitle = "Whole-blood expression associated with mortality",
    x = "log2 Fold Change (Non-Survivor / Survivor)",
    y = "-log10(P-value)",
    color = NULL
  ) +
  theme_pub() +
  ggplot2::theme(legend.position = "top", legend.text = ggplot2::element_text(size = 7.5)) +
  ggplot2::guides(color = ggplot2::guide_legend(nrow = 2))

# Fig 1C: FGSEA Pathway Enrichment Bar/Dot Plot (28-day Non-Survivor vs Survivor)
fg_focus <- dplyr::bind_rows(
  fg_hall_mort %>% dplyr::filter(padj < 0.05) %>% dplyr::slice_head(n = 8),
  fg_reac_mort %>%
    dplyr::filter(grepl("Neutrophil Degranulation|Hypoxia|Glycolysis|Extracellular Matrix|Proteoglycan|Vascular Wall|Hemostasis|Platelet|Fibrin", pathway, ignore.case = TRUE)) %>%
    dplyr::filter(padj < 0.10) %>%
    dplyr::slice_head(n = 7)
) %>%
  dplyr::distinct(pathway, .keep_all = TRUE) %>%
  dplyr::mutate(
    clean_name = gsub("^HALLMARK_", "H: ", pathway),
    clean_name = gsub(" R-HSA-[0-9]+$", "", clean_name),
    clean_name = ifelse(database == "Reactome", paste0("R: ", clean_name), clean_name),
    clean_name = stringr::str_trunc(clean_name, 46),
    direction  = ifelse(NES > 0, "Enriched in 28d Non-Survivors", "Depleted in 28d Non-Survivors")
  ) %>%
  dplyr::arrange(NES) %>%
  dplyr::mutate(clean_name = factor(clean_name, levels = unique(clean_name)))

p1c <- ggplot2::ggplot(fg_focus, ggplot2::aes(x = NES, y = clean_name, fill = direction, size = -log10(pmax(padj, 1e-12)))) +
  ggplot2::geom_segment(ggplot2::aes(x = 0, xend = NES, y = clean_name, yend = clean_name), color = "#666666", linewidth = 0.5) +
  ggplot2::geom_point(shape = 21, color = "#111111", stroke = 0.4) +
  ggplot2::geom_vline(xintercept = 0, linetype = "solid", color = "#222222", linewidth = 0.4) +
  ggplot2::scale_fill_manual(values = c("Enriched in 28d Non-Survivors" = "#B2182B", "Depleted in 28d Non-Survivors" = "#2166AC")) +
  ggplot2::scale_size_continuous(range = c(2.5, 6.0), name = "-log10(FDR)") +
  ggplot2::labs(
    title = "Pathway Dysregulation in 28-Day Non-Survivors",
    subtitle = "FGSEA across Hallmark (H) and Reactome (R) gene sets",
    x = "Normalized Enrichment Score (NES)",
    y = NULL,
    fill = NULL
  ) +
  theme_pub() +
  ggplot2::theme(legend.position = "bottom", legend.box = "vertical")

# Fig 1D: Standardized Trajectory of Locked Signature + Core Vascular Axis Across Clinical States
idx_3grp <- which(pheno65$cohort_group %in% c("Healthy_Control", "NonInfectious_ICU_Control", "Sepsis_28d_Prognostic"))
meta_3grp <- pheno65[idx_3grp, ] %>%
  dplyr::mutate(
    state = dplyr::case_when(
      cohort_group == "Healthy_Control" ~ "Healthy\n(n=42)",
      cohort_group == "NonInfectious_ICU_Control" ~ "Non-Inf ICU\n(n=42)",
      mort_28d == 0 ~ "Sepsis Surv\n(n=365)",
      mort_28d == 1 ~ "Sepsis Non-Surv\n(n=114)"
    ),
    state = factor(state, levels = c("Healthy\n(n=42)", "Non-Inf ICU\n(n=42)", "Sepsis Surv\n(n=365)", "Sepsis Non-Surv\n(n=114)"))
  )

heat_genes <- unique(c(sel_genes, intersect(c("SLC2A1", "LDHA", "ELANE", "MPO", "MMP8", "ARG1", "IL1R2", "F5", "TFPI", "PECAM1"), shared_genes)))
heat_genes <- head(heat_genes, 16)
mat_3grp <- t(scale(t(expr65[heat_genes, idx_3grp, drop = FALSE])))
df_heat <- as.data.frame(mat_3grp) %>%
  tibble::rownames_to_column("gene") %>%
  tidyr::pivot_longer(-gene, names_to = "gsm_id", values_to = "z") %>%
  dplyr::inner_join(meta_3grp[, c("gsm_id", "state")], by = "gsm_id") %>%
  dplyr::group_by(gene, state) %>%
  dplyr::summarise(mean_z = mean(z, na.rm = TRUE), .groups = "drop")

gene_order <- df_heat %>%
  dplyr::filter(state == "Sepsis Non-Surv\n(n=114)") %>%
  dplyr::arrange(mean_z) %>%
  dplyr::pull(gene)
df_heat$gene <- factor(df_heat$gene, levels = gene_order)

p1d <- ggplot2::ggplot(df_heat, ggplot2::aes(x = state, y = gene, fill = mean_z)) +
  ggplot2::geom_tile(color = "white", linewidth = 0.4) +
  ggplot2::geom_text(ggplot2::aes(label = sprintf("%.2f", mean_z)), size = 2.6,
                     color = ifelse(abs(df_heat$mean_z) > 0.65, "white", "#111111"), fontface = "bold") +
  ggplot2::scale_fill_gradient2(low = "#2166AC", mid = "#F7F7F7", high = "#B2182B", midpoint = 0, name = "Mean Z") +
  ggplot2::labs(
    title = "Group-mean blood expression patterns",
    subtitle = "Group-mean standardized expression across GSE65682 states",
    x = NULL,
    y = NULL
  ) +
  theme_pub() +
  ggplot2::theme(axis.text.y = ggplot2::element_text(face = "bold.italic"))

fig1 <- (p1a | p1b) / (p1c | p1d) +
  patchwork::plot_annotation(
    title = "Figure 1. Whole-blood expression and mortality-associated pathway context",
    tag_levels = "A",
    theme = ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", size = 13))
  )

ggplot2::ggsave(file.path(fig_dir, "Fig1_Transcriptomic_Endothelial_Landscape.pdf"), fig1, width = 14, height = 11.5, device = "pdf")
ggplot2::ggsave(file.path(fig_dir, "Fig1_Transcriptomic_Endothelial_Landscape.png"), fig1, width = 14, height = 11.5, dpi = 300)

# --- FIG 2: EIMRS Derivation, Internal Independence & Blind External Validation ---
# Fig 2A: Lollipop plot of Locked LASSO-Cox EIMRS Genes & Weights
sig_plot_df <- sig_table %>%
  dplyr::mutate(
    direction = ifelse(lasso_beta_per_SD > 0, "Mortality Risk (+beta)", "Protective (-beta)"),
    gene = factor(gene, levels = gene[order(lasso_beta_per_SD)])
  )

p2a <- ggplot2::ggplot(sig_plot_df, ggplot2::aes(x = lasso_beta_per_SD, y = gene, color = direction)) +
  ggplot2::geom_vline(xintercept = 0, linetype = "solid", color = "#333333", linewidth = 0.5) +
  ggplot2::geom_segment(ggplot2::aes(x = 0, xend = lasso_beta_per_SD, y = gene, yend = gene), linewidth = 1.0) +
  ggplot2::geom_point(size = 3.4) +
  ggplot2::scale_color_manual(values = c("Mortality Risk (+beta)" = "#B2182B", "Protective (-beta)" = "#2166AC")) +
  ggplot2::labs(
    title = sprintf("Locked %d-Gene EIMRS LASSO-Cox Weights", nrow(sig_table)),
    subtitle = sprintf("GSE65682 (n = 479; full-pipeline bootstrap C = %.3f)", c_opt_corr),
    x = "Penalized LASSO-Cox Coefficient (per 1-SD)",
    y = NULL,
    color = NULL
  ) +
  theme_pub() +
  ggplot2::theme(axis.text.y = ggplot2::element_text(face = "bold.italic"), legend.position = "top")

# Fig 2B: Kaplan-Meier 28-Day Survival by EIMRS Tertiles in GSE65682 (n = 479)
sf_t <- survival::survfit(surv_obj65 ~ eimrs_tertile, data = sepsis_meta65)
lr_t <- survival::survdiff(surv_obj65 ~ eimrs_tertile, data = sepsis_meta65)
p_lr_t <- pchisq(lr_t$chisq, df = length(lr_t$n) - 1, lower.tail = FALSE)

km_df65 <- lapply(seq_along(sf_t$strata), function(i) {
  nm <- sub("^eimrs_tertile=", "", names(sf_t$strata)[i])
  idx_s <- (if (i == 1) 1 else sum(sf_t$strata[1:(i - 1)]) + 1):sum(sf_t$strata[1:i])
  tibble::tibble(
    time = c(0, sf_t$time[idx_s]),
    surv = c(1, sf_t$surv[idx_s]),
    lower = c(1, sf_t$lower[idx_s]),
    upper = c(1, sf_t$upper[idx_s]),
    strata = nm
  )
}) %>% dplyr::bind_rows() %>%
  dplyr::mutate(strata = factor(strata, levels = c("T1 (Low)", "T2 (Intermediate)", "T3 (High)")))

p2b <- ggplot2::ggplot(km_df65, ggplot2::aes(x = time, y = surv, color = strata, fill = strata)) +
  ggplot2::geom_ribbon(ggplot2::aes(ymin = lower, ymax = upper), alpha = 0.14, color = NA) +
  ggplot2::geom_step(linewidth = 1.05) +
  ggplot2::scale_color_manual(values = c("T1 (Low)" = "#2166AC", "T2 (Intermediate)" = "#E08214", "T3 (High)" = "#B2182B")) +
  ggplot2::scale_fill_manual(values = c("T1 (Low)" = "#2166AC", "T2 (Intermediate)" = "#E08214", "T3 (High)" = "#B2182B")) +
  ggplot2::scale_y_continuous(limits = c(0.45, 1.0), labels = scales::percent_format(accuracy = 1)) +
  ggplot2::scale_x_continuous(breaks = seq(0, 28, by = 7)) +
  ggplot2::annotate(
    "text", x = 1.0, y = 0.50, hjust = 0, size = 3.2, fontface = "bold",
    label = sprintf("Log-rank P = %.2e\nHR per 1-SD = %.2f (95%% CI %.2f-%.2f)\nAdj. HR (+MARS/Plt) = %.2f (P = %.2e)",
                    p_lr_t, gse65_cox_summary$HR[1], gse65_cox_summary$CI_Lower[1], gse65_cox_summary$CI_Upper[1],
                    gse65_cox_summary$HR[4], gse65_cox_summary$P_Value[4])
  ) +
  ggplot2::labs(
    title = "28-Day ICU Survival by EIMRS Tertile (GSE65682)",
    subtitle = "GSE65682 development cohort (n = 479, 114 deaths)",
    x = "Days Since ICU Admission",
    y = "28-Day Survival Probability",
    color = "EIMRS Stratum", fill = "EIMRS Stratum"
  ) +
  theme_pub() +
  ggplot2::theme(legend.position = "top")

# Fig 2C: In-Cohort Link to Real Clinical Thrombocytopenia + Blind External GSE95233 Gradient
p2c <- ggplot2::ggplot(meta95_base, ggplot2::aes(x = bio_group, y = eimrs_base_z, fill = bio_group)) +
  ggplot2::geom_violin(alpha = 0.30, color = NA, trim = FALSE) +
  ggplot2::geom_boxplot(width = 0.28, outlier.shape = NA, alpha = 0.85, color = "#111111") +
  ggplot2::geom_jitter(width = 0.12, size = 1.8, alpha = 0.75, color = "#111111") +
  ggplot2::scale_fill_manual(values = c(
    "Healthy\n(n=22)" = "#4DAF4A",
    "Day-1 Shock\nSurvivor (n=34)" = "#2166AC",
    "Day-1 Shock\nNon-Survivor (n=17)" = "#B2182B"
  )) +
  ggplot2::annotate(
    "text", x = 0.65, y = max(meta95_base$eimrs_base_z) + 0.35, hjust = 0, size = 3.1, fontface = "bold",
    label = sprintf("Global Kruskal-Wallis P = %.2e\nDay-1 Surv vs Non-Surv Wilcoxon P = %.4f", kw95_bio_p, wilcox95_p)
  ) +
  ggplot2::labs(
    title = "External Score Distribution (GSE95233 Baseline)",
    subtitle = "22 controls + 51 Day-1 shock samples; pooled cohort scaling",
    x = NULL,
    y = "Locked EIMRS (Standardized Score)"
  ) +
  theme_pub() +
  ggplot2::theme(legend.position = "none")

# Fig 2D: ROC Comparison (GSE65682 Discovery + GSE95233 Blind External Day-1 Validation)
roc_df65 <- tibble::tibble(
  fpr = 1 - roc65$specificities,
  tpr = roc65$sensitivities,
  cohort = sprintf("GSE65682 apparent (n=479): AUC = %.3f (%.3f-%.3f)", auc65_ci[2], auc65_ci[1], auc65_ci[3])
)
roc_df95 <- tibble::tibble(
  fpr = 1 - roc95$specificities,
  tpr = roc95$sensitivities,
  cohort = sprintf("GSE95233 Day-1 Shock (n=51): AUC = %.3f (%.3f-%.3f)", auc95_ci[2], auc95_ci[1], auc95_ci[3])
)
roc_df_oof <- tibble::tibble(
  fpr = 1 - roc_oof$specificities,
  tpr = roc_oof$sensitivities,
  cohort = sprintf("GSE65682 nested OOF (n=479): AUC = %.3f", oof_auc)
)
roc_both <- dplyr::bind_rows(roc_df65, roc_df_oof, roc_df95)

p2d <- ggplot2::ggplot(roc_both, ggplot2::aes(x = fpr, y = tpr, color = cohort)) +
  ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "#888888", linewidth = 0.5) +
  ggplot2::geom_path(linewidth = 1.15) +
  ggplot2::scale_color_manual(values = c("#2166AC", "#439B70", "#B2182B")) +
  ggplot2::annotate(
    "text", x = 0.36, y = 0.18, hjust = 0, size = 3.1, fontface = "bold",
    label = sprintf("External Day-1 OR (per 1-SD):\nUnadj: %.2f (%.2f-%.2f), P = %.4f\nAge+Sex Adj: %.2f (%.2f-%.2f), P = %.4f",
                    or95_uni[1], or95_uni[2], or95_uni[3], s95_uni[4],
                    or95_adj[1], or95_adj[2], or95_adj[3], s95_adj[4])
  ) +
  ggplot2::labs(
    title = "Discrimination in development and external projection",
    subtitle = "Signature-alone external projection; Day-1 cohort-relative scaling",
    x = "1 - Specificity (False Positive Rate)",
    y = "Sensitivity (True Positive Rate)",
    color = NULL
  ) +
  theme_pub() +
  ggplot2::theme(legend.position = "bottom", legend.direction = "vertical",
                 legend.text = ggplot2::element_text(size = 8),
                 legend.background = ggplot2::element_blank())

fig2 <- (p2a | p2b) / (p2c | p2d) +
  patchwork::plot_annotation(
    title = "Figure 2. EIMRS Composition, Survival Associations, and External Projection",
    tag_levels = "A",
    theme = ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", size = 13))
  )

ggplot2::ggsave(file.path(fig_dir, "Fig2_EIMRS_Derivation_and_GSE95233_Validation.pdf"), fig2, width = 14, height = 11.2, device = "pdf")
ggplot2::ggsave(file.path(fig_dir, "Fig2_EIMRS_Derivation_and_GSE95233_Validation.png"), fig2, width = 14, height = 11.2, dpi = 300)

cat("[DONE] 02_transcriptomic_endothelial_signature.R completed successfully.\n")
