#!/usr/bin/env Rscript
# ==============================================================================
# 01_download_and_annotate_geo.R
# Downloads (if absent) and annotates real human ICU whole-blood transcriptomic
# cohorts GSE65682 (MARS Consortium, GPL13667, n = 802; prognostic sepsis n = 479)
# and GSE95233 (Septic Shock Cohort, GPL570, n = 124; Day-1 admission n = 51).
# Maps probes cleanly to harmonized HGNC symbols via EBI A-GEOD-13667 ADF +
# org.Hs.eg.db and hgu133plus2.db, collapses multi-probe genes by maximum IQR,
# and exports reproducible cohort manifests and RDS objects.
# ==============================================================================

suppressPackageStartupMessages({
  library(Biobase)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
  library(hgu133plus2.db)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(tibble)
})

options(timeout = 600)

if (!nzchar(Sys.getenv("CECI_WORK_DIR"))) stop("Set CECI_WORK_DIR to a private working folder.")
proj_dir <- normalizePath(Sys.getenv("CECI_WORK_DIR"), mustWork = TRUE)
raw_dir  <- file.path(proj_dir, "01_data/raw")
proc_dir <- file.path(proj_dir, "01_data/processed")
tab_dir  <- file.path(proj_dir, "03_results/tables")

dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(proc_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)

cat("[1/5] Verifying raw GEO ExpressionSets and GPL13667 ADF annotation...\n")
required_files <- list(
  list(
    path = file.path(raw_dir, "GSE65682_eSet.Rdata"),
    url  = NULL,
    min_bytes = 90000000
  ),
  list(
    path = file.path(raw_dir, "GSE95233_eSet.Rdata"),
    url  = NULL,
    min_bytes = 25000000
  ),
  list(
    path = file.path(raw_dir, "A-GEOD-13667.adf.txt"),
    url  = "https://www.ebi.ac.uk/biostudies/files/A-GEOD-13667/A-GEOD-13667.adf.txt",
    min_bytes = 1500000
  )
)

for (item in required_files) {
  if (!file.exists(item$path) || file.info(item$path)$size < item$min_bytes) {
    if (is.null(item$url)) stop("Missing verified local input: ", item$path)
    cat("  Downloading:", item$url, "-> ", item$path, "\n")
    download.file(item$url, destfile = item$path, mode = "wb", quiet = FALSE)
  } else {
    cat("  Verified existing file:", basename(item$path),
        sprintf("(%.2f MB)\n", file.info(item$path)$size / 1e6))
  }
}

# ------------------------------------------------------------------------------
# 2. Load and curate GSE65682 phenotype metadata
# ------------------------------------------------------------------------------
cat("[2/5] Curating GSE65682 (MARS Consortium) sample phenotypes...\n")
env65 <- new.env()
load(file.path(raw_dir, "GSE65682_eSet.Rdata"), envir = env65)
eset65 <- env65$gset[[1]]
expr65_probe <- Biobase::exprs(eset65)
pd65_raw <- Biobase::pData(eset65)

clean_na <- function(x) {
  x <- trimws(as.character(x))
  ifelse(x == "" | x == "NA" | is.na(x), NA_character_, x)
}

pheno65 <- tibble::tibble(
  gsm_id            = rownames(pd65_raw),
  title             = as.character(pd65_raw$title),
  subject_type_raw  = sub(" \\[.*$", "", as.character(pd65_raw$title)),
  patient_token     = sub("^.*_([^_\\]]+)\\]$", "\\1", as.character(pd65_raw$title)),
  age               = suppressWarnings(as.numeric(clean_na(pd65_raw[["age:ch1"]]))),
  sex               = ifelse(tolower(clean_na(pd65_raw[["gender:ch1"]])) == "male", "Male", "Female"),
  mort_28d          = suppressWarnings(as.integer(clean_na(pd65_raw[["mortality_event_28days:ch1"]]))),
  time_28d_raw      = suppressWarnings(as.numeric(clean_na(pd65_raw[["time_to_event_28days:ch1"]]))),
  mars_endotype     = clean_na(pd65_raw[["endotype_class:ch1"]]),
  mars_subcohort    = clean_na(pd65_raw[["endotype_cohort:ch1"]]),
  pneumonia_dx      = clean_na(pd65_raw[["pneumonia diagnoses:ch1"]]),
  abdo_dx           = clean_na(pd65_raw[["abdominal_sepsis_and_controls:ch1"]]),
  thrombocytopenia  = clean_na(pd65_raw[["thrombocytopenia:ch1"]]),
  diabetes          = clean_na(pd65_raw[["diabetes_mellitus:ch1"]]),
  icuai_status      = clean_na(pd65_raw[["icu_acquired_infection:ch1"]]),
  icuai_paired      = clean_na(pd65_raw[["icu_acquired_infection_paired:ch1"]])
) %>%
  dplyr::mutate(
    # Ensure positive survival time (day 0 deaths coded as 0.5 day for Cox validity)
    time_28d = ifelse(!is.na(time_28d_raw), pmax(0.5, pmin(28.0, time_28d_raw)), NA_real_),
    cohort_group = dplyr::case_when(
      subject_type_raw == "healthy subject" ~ "Healthy_Control",
      abdo_dx == "ctrl_GI" ~ "NonInfectious_ICU_Control",
      !is.na(mort_28d) & !is.na(time_28d) ~ "Sepsis_28d_Prognostic",
      icuai_paired == "Event_ICUAI" ~ "Excluded_Paired_ICUAI_Event",
      TRUE ~ "Other_ICU_Unannotated_Outcome"
    ),
    sepsis_source = dplyr::case_when(
      pneumonia_dx %in% c("cap", "hap") ~ "Pulmonary",
      abdo_dx == "abdo_s" ~ "Abdominal",
      cohort_group == "Sepsis_28d_Prognostic" ~ "Other_Sepsis",
      TRUE ~ NA_character_
    )
  )

cat("  GSE65682 group counts:\n")
print(table(pheno65$cohort_group))
cat("  GSE65682 Sepsis_28d_Prognostic 28-day mortality:\n")
print(table(mort_28d = pheno65$mort_28d[pheno65$cohort_group == "Sepsis_28d_Prognostic"]))

# ------------------------------------------------------------------------------
# 3. Map GPL13667 probes to harmonized HGNC symbols & collapse by max IQR
# ------------------------------------------------------------------------------
cat("[3/5] Mapping GPL13667 (Affymetrix HG-U219) probes via EBI ADF + org.Hs.eg.db...\n")
adf_lines <- readLines(file.path(raw_dir, "A-GEOD-13667.adf.txt"), warn = FALSE)
main_idx <- which(trimws(adf_lines) == "[main]")[1]
adf_df <- read.delim(
  text = paste(adf_lines[(main_idx + 1):length(adf_lines)], collapse = "\n"),
  stringsAsFactors = FALSE,
  check.names = FALSE
)
colnames(adf_df)[1:2] <- c("probe_id", "acc_raw")

acc_long <- adf_df %>%
  dplyr::select(probe_id, acc_raw) %>%
  dplyr::filter(probe_id %in% rownames(expr65_probe), !is.na(acc_raw), acc_raw != "") %>%
  tidyr::separate_rows(acc_raw, sep = ";") %>%
  dplyr::mutate(acc_clean = trimws(sub("\\..*$", "", acc_raw))) %>%
  dplyr::filter(acc_clean != "") %>%
  dplyr::distinct(probe_id, acc_clean)

u_acc <- unique(acc_long$acc_clean)
ref_keys <- intersect(u_acc, AnnotationDbi::keys(org.Hs.eg.db, keytype = "REFSEQ"))
map_ref <- AnnotationDbi::select(
  org.Hs.eg.db, keys = ref_keys, keytype = "REFSEQ", columns = c("ENTREZID", "SYMBOL")
) %>% dplyr::rename(acc_clean = REFSEQ)

rem_acc <- setdiff(u_acc, map_ref$acc_clean[!is.na(map_ref$SYMBOL)])
gb_keys <- intersect(rem_acc, AnnotationDbi::keys(org.Hs.eg.db, keytype = "ACCNUM"))
map_gb <- AnnotationDbi::select(
  org.Hs.eg.db, keys = gb_keys, keytype = "ACCNUM", columns = c("ENTREZID", "SYMBOL")
) %>% dplyr::rename(acc_clean = ACCNUM)

acc_to_gene <- dplyr::bind_rows(map_ref, map_gb) %>%
  dplyr::filter(!is.na(SYMBOL), SYMBOL != "") %>%
  dplyr::distinct(acc_clean, ENTREZID, SYMBOL)

# Keep unambiguous single-gene probes first; if a probe maps to multiple symbols, exclude ambiguous cross-hybridizing probes
probe_gene65 <- acc_long %>%
  dplyr::inner_join(acc_to_gene, by = "acc_clean", relationship = "many-to-many") %>%
  dplyr::distinct(probe_id, SYMBOL, ENTREZID)

ambig_probes65 <- probe_gene65 %>%
  dplyr::count(probe_id) %>%
  dplyr::filter(n > 1) %>%
  dplyr::pull(probe_id)

probe_gene65_clean <- probe_gene65 %>%
  dplyr::filter(!(probe_id %in% ambig_probes65))

# Compute probe IQR across primary discovery/control samples in GSE65682
keep_cols65 <- pheno65$gsm_id[pheno65$cohort_group %in% c("Healthy_Control", "NonInfectious_ICU_Control", "Sepsis_28d_Prognostic")]
iqr65 <- apply(expr65_probe[probe_gene65_clean$probe_id, keep_cols65, drop = FALSE], 1, IQR, na.rm = TRUE)

best_probes65 <- probe_gene65_clean %>%
  dplyr::mutate(probe_iqr = iqr65[probe_id]) %>%
  dplyr::arrange(SYMBOL, dplyr::desc(probe_iqr), probe_id) %>%
  dplyr::group_by(SYMBOL) %>%
  dplyr::slice(1) %>%
  dplyr::ungroup()

expr65_gene <- expr65_probe[best_probes65$probe_id, , drop = FALSE]
rownames(expr65_gene) <- best_probes65$SYMBOL
cat(sprintf("  GSE65682 collapsed to %d unambiguous HGNC genes across %d samples.\n",
            nrow(expr65_gene), ncol(expr65_gene)))

# ------------------------------------------------------------------------------
# 4. Load, curate, and annotate GSE95233 (GPL570, Septic Shock Day-1 + Controls)
# ------------------------------------------------------------------------------
cat("[4/5] Curating and annotating GSE95233 (GPL570) baseline cohort...\n")
env95 <- new.env()
load(file.path(raw_dir, "GSE95233_eSet.Rdata"), envir = env95)
eset95 <- env95$gset[[1]]
expr95_probe <- Biobase::exprs(eset95)
pd95_raw <- Biobase::pData(eset95)

pheno95 <- tibble::tibble(
  gsm_id          = rownames(pd95_raw),
  title           = as.character(pd95_raw$title),
  source_name     = as.character(pd95_raw$source_name_ch1),
  age             = suppressWarnings(as.numeric(clean_na(pd95_raw[["age:ch1"]]))),
  gender_raw      = clean_na(pd95_raw[["gender:ch1"]]),
  # French annotation 'H' = Homme (Male), 'M' = Male, 'F' = Female
  sex             = ifelse(toupper(clean_na(pd95_raw[["gender:ch1"]])) %in% c("M", "H"), "Male", "Female"),
  time_point      = clean_na(pd95_raw[["time point:ch1"]]),
  survival_raw    = clean_na(pd95_raw[["survival:ch1"]])
) %>%
  dplyr::mutate(
    cohort_group = dplyr::case_when(
      is.na(time_point) & grepl("^Control", source_name, ignore.case = TRUE) ~ "Healthy_Control",
      time_point == "D01" & survival_raw %in% c("Survivor", "Non Survivor") ~ "Septic_Shock_Day1",
      time_point %in% c("D02", "D03") ~ "Excluded_Serial_D02_D03",
      TRUE ~ "Other"
    ),
    mort_event = dplyr::case_when(
      survival_raw == "Non Survivor" ~ 1L,
      survival_raw == "Survivor" ~ 0L,
      TRUE ~ NA_integer_
    )
  )

cat("  GSE95233 group counts:\n")
print(table(pheno95$cohort_group))
cat("  GSE95233 Day-1 Septic Shock mortality:\n")
print(table(mort_event = pheno95$mort_event[pheno95$cohort_group == "Septic_Shock_Day1"]))

# Map GPL570 via hgu133plus2.db ENTREZID -> org.Hs.eg.db SYMBOL for 100% harmonized nomenclature
map95_entrez <- AnnotationDbi::select(
  hgu133plus2.db,
  keys = rownames(expr95_probe),
  keytype = "PROBEID",
  columns = "ENTREZID"
) %>%
  dplyr::filter(!is.na(ENTREZID)) %>%
  dplyr::distinct(PROBEID, ENTREZID)

entrez_sym <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = unique(map95_entrez$ENTREZID),
  keytype = "ENTREZID",
  columns = "SYMBOL"
) %>%
  dplyr::filter(!is.na(SYMBOL), SYMBOL != "") %>%
  dplyr::distinct(ENTREZID, SYMBOL)

probe_gene95 <- map95_entrez %>%
  dplyr::inner_join(entrez_sym, by = "ENTREZID") %>%
  dplyr::distinct(PROBEID, SYMBOL, ENTREZID)

ambig_probes95 <- probe_gene95 %>%
  dplyr::count(PROBEID) %>%
  dplyr::filter(n > 1) %>%
  dplyr::pull(PROBEID)

probe_gene95_clean <- probe_gene95 %>%
  dplyr::filter(!(PROBEID %in% ambig_probes95))

# Compute probe IQR strictly on baseline samples (Healthy_Control + Septic_Shock_Day1) - zero D02/D03 leakage
base_cols95 <- pheno95$gsm_id[pheno95$cohort_group %in% c("Healthy_Control", "Septic_Shock_Day1")]
iqr95 <- apply(expr95_probe[probe_gene95_clean$PROBEID, base_cols95, drop = FALSE], 1, IQR, na.rm = TRUE)

best_probes95 <- probe_gene95_clean %>%
  dplyr::mutate(probe_iqr = iqr95[PROBEID]) %>%
  dplyr::arrange(SYMBOL, dplyr::desc(probe_iqr), PROBEID) %>%
  dplyr::group_by(SYMBOL) %>%
  dplyr::slice(1) %>%
  dplyr::ungroup()

expr95_gene <- expr95_probe[best_probes95$PROBEID, , drop = FALSE]
rownames(expr95_gene) <- best_probes95$SYMBOL
cat(sprintf("  GSE95233 collapsed to %d unambiguous HGNC genes across %d samples.\n",
            nrow(expr95_gene), ncol(expr95_gene)))

# ------------------------------------------------------------------------------
# 5. Export harmonized manifests, summary table, and RDS bundle
# ------------------------------------------------------------------------------
cat("[5/5] Exporting processed manifests and RDS bundle...\n")
shared_genes <- sort(intersect(rownames(expr65_gene), rownames(expr95_gene)))
cat(sprintf("  Harmonized 1:1 shared HGNC gene universe: %d genes.\n", length(shared_genes)))

readr::write_csv(pheno65, file.path(proc_dir, "gse65682_pheno_clean.csv"))
readr::write_csv(pheno95, file.path(proc_dir, "gse95233_pheno_clean.csv"))

cohort_summary <- tibble::tibble(
  Cohort_ID = c("GSE65682 (MARS Discovery)", "GSE65682 (MARS Controls)", "GSE95233 (External Validation)", "GSE95233 (External Controls)"),
  Platform = c("GPL13667 (Affymetrix HG-U219)", "GPL13667 (Affymetrix HG-U219)", "GPL570 (Affymetrix HG-U133 Plus 2.0)", "GPL570 (Affymetrix HG-U133 Plus 2.0)"),
  Population = c("Adult ICU Sepsis (28-day survival follow-up)", "Healthy Volunteers (n=42) + Non-Infectious GI ICU Controls (n=42)", "Adult ICU Septic Shock (Day-1 baseline admission only)", "Healthy Volunteers"),
  Sample_Size_N = c(
    sum(pheno65$cohort_group == "Sepsis_28d_Prognostic"),
    sum(pheno65$cohort_group %in% c("Healthy_Control", "NonInfectious_ICU_Control")),
    sum(pheno95$cohort_group == "Septic_Shock_Day1"),
    sum(pheno95$cohort_group == "Healthy_Control")
  ),
  Mortality_Events_n_pct = c(
    sprintf("%d (%.1f%%)", sum(pheno65$mort_28d == 1, na.rm = TRUE), 100 * mean(pheno65$mort_28d[pheno65$cohort_group == "Sepsis_28d_Prognostic"] == 1)),
    "0 (0.0%)",
    sprintf("%d (%.1f%%)", sum(pheno95$mort_event[pheno95$cohort_group == "Septic_Shock_Day1"] == 1), 100 * mean(pheno95$mort_event[pheno95$cohort_group == "Septic_Shock_Day1"] == 1)),
    "0 (0.0%)"
  ),
  Mapped_Unambiguous_Genes = c(nrow(expr65_gene), nrow(expr65_gene), nrow(expr95_gene), nrow(expr95_gene)),
  Shared_Gene_Universe = length(shared_genes)
)
readr::write_csv(cohort_summary, file.path(tab_dir, "Table_S0_Transcriptomic_Cohort_Summary.csv"))
print(cohort_summary)

saveRDS(
  list(
    expr65_gene   = expr65_gene,
    pheno65       = pheno65,
    probe_map65   = best_probes65,
    expr95_gene   = expr95_gene,
    pheno95       = pheno95,
    probe_map95   = best_probes95,
    shared_genes  = shared_genes
  ),
  file = file.path(proc_dir, "geo_annotated_cohorts.rds")
)
cat("[DONE] 01_download_and_annotate_geo.R completed successfully.\n")
