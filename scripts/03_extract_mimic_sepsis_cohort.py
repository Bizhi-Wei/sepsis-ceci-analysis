#!/usr/bin/env python3
"""
03_extract_mimic_sepsis_cohort.py
=================================
Read-only extraction of the adult ICU Sepsis cohort from MIMIC-IV v3.1 (the configured MIMIC_DIR)
using in-memory DuckDB. Zero files are created or modified inside the configured MIMIC_DIR.

Cohort Architecture (24-Hour Landmark Design to Eliminate Immortal-Time Bias):
  1. All ICU stays in MIMIC-IV v3.1 (icustays.csv)
  2. Adult patients aged >= 18 years at ICU admission
  3. First ICU stay per unique subject_id (eliminates repeated-measures clustering)
  4. ICU length of stay >= 24 hours AND alive > 24 hours from ICU intime (24h landmark)
  5. Clinical Sepsis present at ICU admission: explicit Sepsis ICD-9/10 diagnosis
     WITH active microbiological culture workup within [-48h, +24h] of ICU intime
     (prevents late ICU-acquired nosocomial sepsis misclassification at Day 1)
  6. Primary Analytic Cohort: complete Day-1 [-6h, +24h] endothelial-coagulopathy
     biomarker panel (Lactate, Serum Albumin, Platelet Count, INR).
     Full eligible cohort (Step 5) is also saved for missingness sensitivity & MICE/imputation checks.
"""

import time
from pathlib import Path
import os
import duckdb
import numpy as np
import pandas as pd

PROJECT_DIR = Path(os.environ["CECI_WORK_DIR"]).resolve(strict=True)
MIMIC_DIR = Path(os.environ["MIMIC_DIR"]).resolve(strict=True)
PROC_DIR = PROJECT_DIR / "01_data" / "processed"
TAB_DIR = PROJECT_DIR / "03_results" / "tables"
PROC_DIR.mkdir(parents=True, exist_ok=True)
TAB_DIR.mkdir(parents=True, exist_ok=True)


def main():
    t0 = time.time()
    print("[1/6] Initializing in-memory DuckDB connection (100% READ-ONLY on the configured MIMIC_DIR)...")
    con = duckdb.connect(":memory:")
    con.execute("SET threads=8; SET memory_limit='6GB'")

    # Step 1-5: Build base ICU stays, STROBE attrition, and eligible ICU-admission Sepsis cohort
    print("[2/6] Screening ICU stays, patient demographics, sepsis diagnoses, and [-48h, +24h] cultures...")
    base_sql = f"""
    CREATE TEMP TABLE all_stays AS
    SELECT
        i.subject_id,
        i.hadm_id,
        i.stay_id,
        i.first_careunit,
        CAST(i.intime AS TIMESTAMP) AS intime,
        CAST(i.outtime AS TIMESTAMP) AS outtime,
        CAST(i.los AS DOUBLE) AS icu_los_days,
        p.gender,
        CAST(p.anchor_age AS INT) + EXTRACT(YEAR FROM CAST(i.intime AS TIMESTAMP)) - CAST(p.anchor_year AS INT) AS age,
        CAST(p.dod AS TIMESTAMP) AS dod,
        CAST(a.admittime AS TIMESTAMP) AS admittime,
        CAST(a.dischtime AS TIMESTAMP) AS dischtime,
        CAST(a.deathtime AS TIMESTAMP) AS deathtime,
        a.admission_type,
        a.race,
        CAST(a.hospital_expire_flag AS INT) AS hospital_expire_flag,
        ROW_NUMBER() OVER (PARTITION BY i.subject_id ORDER BY CAST(i.intime AS TIMESTAMP) ASC, i.stay_id ASC) AS stay_rn
    FROM read_csv_auto('{MIMIC_DIR.as_posix()}/icu/icustays.csv') i
    JOIN read_csv_auto('{MIMIC_DIR.as_posix()}/hosp/patients.csv') p ON i.subject_id = p.subject_id
    JOIN read_csv_auto('{MIMIC_DIR.as_posix()}/hosp/admissions.csv') a ON i.hadm_id = a.hadm_id;
    """
    con.execute(base_sql)

    # Comorbidities & Sepsis diagnosis extraction from diagnoses_icd
    dx_sql = f"""
    CREATE TEMP TABLE hadm_dx AS
    SELECT
        hadm_id,
        MAX(CASE
            WHEN (icd_version = 9 AND (icd_code IN ('99591','99592','78552') OR icd_code LIKE '038%'))
              OR (icd_version = 10 AND (icd_code LIKE 'A40%' OR icd_code LIKE 'A41%' OR icd_code IN ('R6520','R6521')))
            THEN 1 ELSE 0 END) AS sepsis_icd,
        MAX(CASE
            WHEN (icd_version = 9 AND icd_code = '78552')
              OR (icd_version = 10 AND (icd_code = 'R6521' OR icd_code LIKE 'T8112%'))
            THEN 1 ELSE 0 END) AS septic_shock_icd,
        MAX(CASE
            WHEN (icd_version = 9 AND (icd_code LIKE '428%' OR icd_code IN ('39891','40201','40211','40291','40401','40403','40411','40413','40491','40493')))
              OR (icd_version = 10 AND (icd_code LIKE 'I50%' OR icd_code IN ('I099','I110','I130','I132','I255','I420','I425','I426','I427','I428','I429','P290')))
            THEN 1 ELSE 0 END) AS chf,
        MAX(CASE
            WHEN (icd_version = 9 AND (icd_code LIKE '585%' OR icd_code LIKE '582%' OR icd_code LIKE '583%' OR icd_code IN ('586','5880','V420') OR icd_code LIKE 'V451%' OR icd_code LIKE 'V56%' OR icd_code LIKE '403%'))
              OR (icd_version = 10 AND (icd_code LIKE 'N18%' OR icd_code IN ('I120','I131','N19','N250','Z490','Z491','Z492','Z940','Z992')))
            THEN 1 ELSE 0 END) AS ckd,
        MAX(CASE
            WHEN (icd_version = 9 AND (icd_code LIKE '571%' OR icd_code IN ('07022','07023','07032','07033','07044','07054','0706','0709','570','5733','5734','5738','5739','V427','4560','4561','4562','5722','5723','5724','5728')))
              OR (icd_version = 10 AND (icd_code LIKE 'B18%' OR icd_code LIKE 'K70%' OR icd_code LIKE 'K73%' OR icd_code LIKE 'K74%' OR icd_code IN ('K713','K714','K715','K717','K760','K762','K763','K764','K768','K769','Z944','I850','I859','I864','I982','K721','K729','K765','K766','K767')))
            THEN 1 ELSE 0 END) AS liver_disease,
        MAX(CASE
            WHEN (icd_version = 9 AND icd_code IN ('4560','4561','4562','5722','5723','5724','5728','5712','5715','5716'))
              OR (icd_version = 10 AND (icd_code IN ('I850','I859','I864','I982','K704','K711','K721','K729','K765','K766','K767') OR icd_code LIKE 'K74%'))
            THEN 1 ELSE 0 END) AS cirrhosis_portal_htn,
        MAX(CASE
            WHEN (icd_version = 9 AND icd_code LIKE '250%')
              OR (icd_version = 10 AND (icd_code LIKE 'E10%' OR icd_code LIKE 'E11%' OR icd_code LIKE 'E13%'))
            THEN 1 ELSE 0 END) AS diabetes,
        MAX(CASE
            WHEN (icd_version = 9 AND (icd_code BETWEEN '140' AND '209' OR icd_code LIKE '14%' OR icd_code LIKE '15%' OR icd_code LIKE '16%' OR icd_code LIKE '17%' OR icd_code LIKE '18%' OR icd_code LIKE '19%' OR icd_code LIKE '20%'))
              OR (icd_version = 10 AND (icd_code LIKE 'C0%' OR icd_code LIKE 'C1%' OR icd_code LIKE 'C2%' OR icd_code LIKE 'C3%' OR icd_code LIKE 'C4%' OR icd_code LIKE 'C5%' OR icd_code LIKE 'C6%' OR icd_code LIKE 'C7%' OR icd_code LIKE 'C8%' OR icd_code LIKE 'C9%'))
            THEN 1 ELSE 0 END) AS malignancy,
        MAX(CASE
            WHEN (icd_version = 9 AND (icd_code LIKE '490%' OR icd_code LIKE '491%' OR icd_code LIKE '492%' OR icd_code LIKE '493%' OR icd_code LIKE '494%' OR icd_code LIKE '495%' OR icd_code LIKE '496%'))
              OR (icd_version = 10 AND (icd_code LIKE 'J40%' OR icd_code LIKE 'J41%' OR icd_code LIKE 'J42%' OR icd_code LIKE 'J43%' OR icd_code LIKE 'J44%' OR icd_code LIKE 'J45%' OR icd_code LIKE 'J47%'))
            THEN 1 ELSE 0 END) AS copd,
        MAX(CASE
            WHEN (icd_version = 9 AND (icd_code LIKE '410%' OR icd_code LIKE '412%'))
              OR (icd_version = 10 AND (icd_code LIKE 'I21%' OR icd_code LIKE 'I22%' OR icd_code = 'I252'))
            THEN 1 ELSE 0 END) AS cad_mi,
        MAX(CASE
            WHEN (icd_version = 9 AND (icd_code LIKE '430%' OR icd_code LIKE '431%' OR icd_code LIKE '432%' OR icd_code LIKE '433%' OR icd_code LIKE '434%' OR icd_code LIKE '435%' OR icd_code LIKE '436%' OR icd_code LIKE '437%' OR icd_code LIKE '438%'))
              OR (icd_version = 10 AND (icd_code LIKE 'I60%' OR icd_code LIKE 'I61%' OR icd_code LIKE 'I62%' OR icd_code LIKE 'I63%' OR icd_code LIKE 'I64%' OR icd_code LIKE 'I65%' OR icd_code LIKE 'I66%' OR icd_code LIKE 'I67%' OR icd_code LIKE 'I69%' OR icd_code LIKE 'G45%'))
            THEN 1 ELSE 0 END) AS cerebrovascular,
        MAX(CASE
            WHEN (icd_version = 9 AND (icd_code LIKE '480%' OR icd_code LIKE '481%' OR icd_code LIKE '482%' OR icd_code LIKE '483%' OR icd_code LIKE '485%' OR icd_code LIKE '486%'))
              OR (icd_version = 10 AND (icd_code LIKE 'J12%' OR icd_code LIKE 'J13%' OR icd_code LIKE 'J14%' OR icd_code LIKE 'J15%' OR icd_code LIKE 'J16%' OR icd_code LIKE 'J17%' OR icd_code LIKE 'J18%'))
            THEN 1 ELSE 0 END) AS source_pulmonary,
        MAX(CASE
            WHEN (icd_version = 9 AND (icd_code LIKE '540%' OR icd_code LIKE '562%' OR icd_code LIKE '567%' OR icd_code LIKE '575%' OR icd_code LIKE '576%' OR icd_code LIKE '008%'))
              OR (icd_version = 10 AND (icd_code LIKE 'K35%' OR icd_code LIKE 'K57%' OR icd_code LIKE 'K65%' OR icd_code LIKE 'K81%' OR icd_code LIKE 'K830%' OR icd_code LIKE 'A04%'))
            THEN 1 ELSE 0 END) AS source_abdominal,
        MAX(CASE
            WHEN (icd_version = 9 AND (icd_code LIKE '590%' OR icd_code = '5990'))
              OR (icd_version = 10 AND (icd_code LIKE 'N10%' OR icd_code LIKE 'N12%' OR icd_code = 'N390'))
            THEN 1 ELSE 0 END) AS source_urinary
    FROM read_csv_auto('{MIMIC_DIR.as_posix()}/hosp/diagnoses_icd.csv')
    GROUP BY hadm_id;
    """
    con.execute(dx_sql)

    # Active microbiological culture around ICU admission [-48h, +24h]
    micro_sql = f"""
    CREATE TEMP TABLE early_micro AS
    SELECT
        st.stay_id,
        1 AS early_culture_flag,
        MAX(CASE WHEN m.org_name IS NOT NULL AND TRIM(CAST(m.org_name AS VARCHAR)) <> '' THEN 1 ELSE 0 END) AS culture_positive_flag
    FROM all_stays st
    JOIN hadm_dx dx ON st.hadm_id = dx.hadm_id AND dx.sepsis_icd = 1
    JOIN read_csv_auto('{MIMIC_DIR.as_posix()}/hosp/microbiologyevents.csv') m ON st.hadm_id = m.hadm_id
    WHERE st.age >= 18
      AND st.stay_rn = 1
      AND st.icu_los_days >= 1.0
      AND (st.deathtime IS NULL OR st.deathtime > st.intime + INTERVAL 24 HOUR)
      AND COALESCE(TRY_CAST(m.charttime AS TIMESTAMP), TRY_CAST(m.chartdate AS TIMESTAMP))
          BETWEEN st.intime - INTERVAL 48 HOUR AND st.intime + INTERVAL 24 HOUR
    GROUP BY st.stay_id;
    """
    con.execute(micro_sql)

    # Compute STROBE attrition up to Step 5
    strobe_pre = con.execute("""
    SELECT
        COUNT(*) AS step1_all_icu_stays,
        SUM(CASE WHEN age >= 18 THEN 1 ELSE 0 END) AS step2_adult_stays,
        SUM(CASE WHEN age >= 18 AND stay_rn = 1 THEN 1 ELSE 0 END) AS step3_first_icu_stay,
        SUM(CASE WHEN age >= 18 AND stay_rn = 1 AND icu_los_days >= 1.0
                  AND (deathtime IS NULL OR deathtime > intime + INTERVAL 24 HOUR)
                 THEN 1 ELSE 0 END) AS step4_los_ge_24h_landmark,
        SUM(CASE WHEN age >= 18 AND stay_rn = 1 AND icu_los_days >= 1.0
                  AND (deathtime IS NULL OR deathtime > intime + INTERVAL 24 HOUR)
                  AND COALESCE(dx.sepsis_icd, 0) = 1
                 THEN 1 ELSE 0 END) AS step5a_sepsis_icd_any,
        SUM(CASE WHEN age >= 18 AND stay_rn = 1 AND icu_los_days >= 1.0
                  AND (deathtime IS NULL OR deathtime > intime + INTERVAL 24 HOUR)
                  AND COALESCE(dx.sepsis_icd, 0) = 1
                  AND COALESCE(em.early_culture_flag, 0) = 1
                 THEN 1 ELSE 0 END) AS step5_sepsis_admission_verified
    FROM all_stays st
    LEFT JOIN hadm_dx dx ON st.hadm_id = dx.hadm_id
    LEFT JOIN early_micro em ON st.stay_id = em.stay_id
    """).df().iloc[0].to_dict()
    print("  STROBE pre-biomarker counts:", strobe_pre)

    # Create eligible_sepsis table (Step 5: N = 7,195)
    con.execute("""
    CREATE TEMP TABLE eligible_sepsis AS
    SELECT
        st.*,
        COALESCE(dx.septic_shock_icd, 0) AS septic_shock_icd,
        COALESCE(dx.chf, 0) AS chf,
        COALESCE(dx.ckd, 0) AS ckd,
        COALESCE(dx.liver_disease, 0) AS liver_disease,
        COALESCE(dx.cirrhosis_portal_htn, 0) AS cirrhosis_portal_htn,
        COALESCE(dx.diabetes, 0) AS diabetes,
        COALESCE(dx.malignancy, 0) AS malignancy,
        COALESCE(dx.copd, 0) AS copd,
        COALESCE(dx.cad_mi, 0) AS cad_mi,
        COALESCE(dx.cerebrovascular, 0) AS cerebrovascular,
        COALESCE(dx.source_pulmonary, 0) AS source_pulmonary,
        COALESCE(dx.source_abdominal, 0) AS source_abdominal,
        COALESCE(dx.source_urinary, 0) AS source_urinary,
        COALESCE(em.culture_positive_flag, 0) AS culture_positive_flag
    FROM all_stays st
    JOIN hadm_dx dx ON st.hadm_id = dx.hadm_id AND dx.sepsis_icd = 1
    JOIN early_micro em ON st.stay_id = em.stay_id AND em.early_culture_flag = 1
    WHERE st.age >= 18
      AND st.stay_rn = 1
      AND st.icu_los_days >= 1.0
      AND (st.deathtime IS NULL OR st.deathtime > st.intime + INTERVAL 24 HOUR);
    """)

    # Step 3/6: Extract ICU organ support (vasopressors, IMV, RRT) in [intime, intime + 24h]
    print(f"[3/6] Extracting Day-1 vasopressors, invasive mechanical ventilation, and RRT ({time.time()-t0:.1f}s)...")
    con.execute(f"""
    CREATE TEMP TABLE pressors_24h AS
    SELECT
        es.stay_id,
        1 AS vasopressor_24h,
        MAX(CASE WHEN TRY_CAST(ie.itemid AS BIGINT) = 221906 THEN 1 ELSE 0 END) AS norepi_24h,
        MAX(CASE WHEN TRY_CAST(ie.itemid AS BIGINT) = 222315 THEN 1 ELSE 0 END) AS vasopressin_24h
    FROM eligible_sepsis es
    JOIN read_csv('{MIMIC_DIR.as_posix()}/icu/inputevents.csv', header=true, all_varchar=true, sample_size=20480) ie
      ON es.stay_id = TRY_CAST(ie.stay_id AS BIGINT)
    WHERE TRY_CAST(ie.itemid AS BIGINT) IN (221906, 222315, 221289, 229617, 221662, 221749, 229630, 229631, 229632, 229709, 229764)
      AND (TRY_CAST(ie.amount AS DOUBLE) > 0 OR TRY_CAST(ie.rate AS DOUBLE) > 0)
      AND TRY_CAST(ie.starttime AS TIMESTAMP) < es.intime + INTERVAL 24 HOUR
      AND TRY_CAST(ie.endtime AS TIMESTAMP) > es.intime
    GROUP BY es.stay_id;
    """)

    con.execute(f"""
    CREATE TEMP TABLE proc_24h AS
    SELECT
        es.stay_id,
        MAX(CASE WHEN TRY_CAST(pe.itemid AS BIGINT) IN (225792, 224385) THEN 1 ELSE 0 END) AS imv_24h,
        MAX(CASE WHEN TRY_CAST(pe.itemid AS BIGINT) IN (225802, 225803, 225805, 225809, 225955, 225441) THEN 1 ELSE 0 END) AS rrt_24h
    FROM eligible_sepsis es
    JOIN read_csv('{MIMIC_DIR.as_posix()}/icu/procedureevents.csv', header=true, all_varchar=true, sample_size=20480) pe
      ON es.stay_id = TRY_CAST(pe.stay_id AS BIGINT)
    WHERE TRY_CAST(pe.itemid AS BIGINT) IN (225792, 224385, 225802, 225803, 225805, 225809, 225955, 225441)
      AND TRY_CAST(pe.starttime AS TIMESTAMP) < es.intime + INTERVAL 24 HOUR
      AND TRY_CAST(pe.endtime AS TIMESTAMP) > es.intime
    GROUP BY es.stay_id;
    """)

    # Step 4/6: Single-pass extraction of Day-1 labs from hosp/labevents.csv
    print(f"[4/6] Scanning hosp/labevents.csv in single pass for Day-1 [-6h, +24h] & [-24h, +48h] biomarkers ({time.time()-t0:.1f}s)...")
    con.execute("""
    CREATE TEMP TABLE cand_hadm AS SELECT DISTINCT hadm_id FROM eligible_sepsis;
    """)
    con.execute(f"""
    CREATE TEMP TABLE raw_labs AS
    SELECT
        TRY_CAST(l.hadm_id AS BIGINT) AS hadm_id,
        TRY_CAST(l.itemid AS BIGINT) AS itemid,
        TRY_CAST(l.charttime AS TIMESTAMP) AS charttime,
        TRY_CAST(l.valuenum AS DOUBLE) AS val
    FROM read_csv('{MIMIC_DIR.as_posix()}/hosp/labevents.csv', header=true, all_varchar=true, sample_size=20480) l
    SEMI JOIN cand_hadm ch ON TRY_CAST(l.hadm_id AS BIGINT) = ch.hadm_id
    WHERE TRY_CAST(l.itemid AS BIGINT) IN (
        50813, 52442, 53154,        -- Lactate
        50862, 53085, 53138, 52022, -- Albumin
        51265, 53189,               -- Platelets
        51237, 51675,               -- INR
        51300, 51301, 51755, 51756, -- WBC
        50912, 52546, 52024,        -- Creatinine
        50885, 53089,               -- Total Bilirubin
        51006, 52647,               -- BUN
        50882, 50803,               -- Bicarbonate
        51222, 50811, 51640         -- Hemoglobin
    )
    AND TRY_CAST(l.valuenum AS DOUBLE) IS NOT NULL;
    """)

    con.execute("""
    CREATE TEMP TABLE labs_agg AS
    SELECT
        es.stay_id,
        -- Primary [-6h, +24h] window around ICU intime
        MAX(CASE WHEN rl.itemid IN (50813,52442,53154) AND rl.val BETWEEN 0.1 AND 30.0
                  AND rl.charttime BETWEEN es.intime - INTERVAL 6 HOUR AND es.intime + INTERVAL 24 HOUR
                 THEN rl.val END) AS lactate_24h,
        MIN(CASE WHEN rl.itemid IN (50862,53085,53138,52022) AND rl.val BETWEEN 0.5 AND 6.5
                  AND rl.charttime BETWEEN es.intime - INTERVAL 6 HOUR AND es.intime + INTERVAL 24 HOUR
                 THEN rl.val END) AS albumin_24h,
        MIN(CASE WHEN rl.itemid IN (51265,53189) AND rl.val BETWEEN 1.0 AND 1500.0
                  AND rl.charttime BETWEEN es.intime - INTERVAL 6 HOUR AND es.intime + INTERVAL 24 HOUR
                 THEN rl.val END) AS platelet_24h,
        MAX(CASE WHEN rl.itemid IN (51237,51675) AND rl.val BETWEEN 0.5 AND 25.0
                  AND rl.charttime BETWEEN es.intime - INTERVAL 6 HOUR AND es.intime + INTERVAL 24 HOUR
                 THEN rl.val END) AS inr_24h,
        MAX(CASE WHEN rl.itemid IN (51300,51301,51755,51756) AND rl.val BETWEEN 0.1 AND 200.0
                  AND rl.charttime BETWEEN es.intime - INTERVAL 6 HOUR AND es.intime + INTERVAL 24 HOUR
                 THEN rl.val END) AS wbc_24h,
        MAX(CASE WHEN rl.itemid IN (50912,52546,52024) AND rl.val BETWEEN 0.1 AND 25.0
                  AND rl.charttime BETWEEN es.intime - INTERVAL 6 HOUR AND es.intime + INTERVAL 24 HOUR
                 THEN rl.val END) AS creatinine_24h,
        MAX(CASE WHEN rl.itemid IN (50885,53089) AND rl.val BETWEEN 0.05 AND 60.0
                  AND rl.charttime BETWEEN es.intime - INTERVAL 6 HOUR AND es.intime + INTERVAL 24 HOUR
                 THEN rl.val END) AS bilirubin_24h,
        MAX(CASE WHEN rl.itemid IN (51006,52647) AND rl.val BETWEEN 1.0 AND 250.0
                  AND rl.charttime BETWEEN es.intime - INTERVAL 6 HOUR AND es.intime + INTERVAL 24 HOUR
                 THEN rl.val END) AS bun_24h,
        MIN(CASE WHEN rl.itemid IN (50882,50803) AND rl.val BETWEEN 2.0 AND 55.0
                  AND rl.charttime BETWEEN es.intime - INTERVAL 6 HOUR AND es.intime + INTERVAL 24 HOUR
                 THEN rl.val END) AS bicarbonate_24h,
        MIN(CASE WHEN rl.itemid IN (51222,50811,51640) AND rl.val BETWEEN 2.0 AND 25.0
                  AND rl.charttime BETWEEN es.intime - INTERVAL 6 HOUR AND es.intime + INTERVAL 24 HOUR
                 THEN rl.val END) AS hemoglobin_24h,
        -- Extended [-24h, +24h] pre-ICU ED/ward admission window for sensitivity check
        MAX(CASE WHEN rl.itemid IN (50813,52442,53154) AND rl.val BETWEEN 0.1 AND 30.0
                  AND rl.charttime BETWEEN es.intime - INTERVAL 24 HOUR AND es.intime + INTERVAL 24 HOUR
                 THEN rl.val END) AS lactate_ext24h,
        MIN(CASE WHEN rl.itemid IN (50862,53085,53138,52022) AND rl.val BETWEEN 0.5 AND 6.5
                  AND rl.charttime BETWEEN es.intime - INTERVAL 24 HOUR AND es.intime + INTERVAL 24 HOUR
                 THEN rl.val END) AS albumin_ext24h,
        MIN(CASE WHEN rl.itemid IN (51265,53189) AND rl.val BETWEEN 1.0 AND 1500.0
                  AND rl.charttime BETWEEN es.intime - INTERVAL 24 HOUR AND es.intime + INTERVAL 24 HOUR
                 THEN rl.val END) AS platelet_ext24h,
        MAX(CASE WHEN rl.itemid IN (51237,51675) AND rl.val BETWEEN 0.5 AND 25.0
                  AND rl.charttime BETWEEN es.intime - INTERVAL 24 HOUR AND es.intime + INTERVAL 24 HOUR
                 THEN rl.val END) AS inr_ext24h
    FROM eligible_sepsis es
    LEFT JOIN raw_labs rl ON es.hadm_id = rl.hadm_id
    GROUP BY es.stay_id;
    """)

    # Step 5/6: Single-pass extraction of Day-1 vital signs from icu/chartevents.csv
    print(f"[5/6] Scanning icu/chartevents.csv in single pass for Day-1 MAP, HR, RR, SpO2, Temp ({time.time()-t0:.1f}s)...")
    con.execute("""
    CREATE TEMP TABLE cand_stays AS SELECT DISTINCT stay_id FROM eligible_sepsis;
    """)
    con.execute(f"""
    CREATE TEMP TABLE raw_vitals AS
    SELECT
        TRY_CAST(c.stay_id AS BIGINT) AS stay_id,
        TRY_CAST(c.itemid AS BIGINT) AS itemid,
        TRY_CAST(c.charttime AS TIMESTAMP) AS charttime,
        TRY_CAST(c.valuenum AS DOUBLE) AS val
    FROM read_csv('{MIMIC_DIR.as_posix()}/icu/chartevents.csv', header=true, all_varchar=true, sample_size=20480) c
    SEMI JOIN cand_stays cs ON TRY_CAST(c.stay_id AS BIGINT) = cs.stay_id
    WHERE TRY_CAST(c.itemid AS BIGINT) IN (220052, 220181, 225312, 220045, 220210, 224690, 220277, 223761, 223762)
      AND TRY_CAST(c.valuenum AS DOUBLE) IS NOT NULL;
    """)

    con.execute("""
    CREATE TEMP TABLE vitals_agg AS
    SELECT
        es.stay_id,
        MIN(CASE WHEN rv.itemid IN (220052,220181,225312) AND rv.val BETWEEN 20.0 AND 200.0
                  AND rv.charttime BETWEEN es.intime AND es.intime + INTERVAL 24 HOUR
                 THEN rv.val END) AS map_min_24h,
        AVG(CASE WHEN rv.itemid IN (220052,220181,225312) AND rv.val BETWEEN 20.0 AND 200.0
                  AND rv.charttime BETWEEN es.intime AND es.intime + INTERVAL 24 HOUR
                 THEN rv.val END) AS map_mean_24h,
        MAX(CASE WHEN rv.itemid = 220045 AND rv.val BETWEEN 20.0 AND 250.0
                  AND rv.charttime BETWEEN es.intime AND es.intime + INTERVAL 24 HOUR
                 THEN rv.val END) AS hr_max_24h,
        MAX(CASE WHEN rv.itemid IN (220210,224690) AND rv.val BETWEEN 4.0 AND 70.0
                  AND rv.charttime BETWEEN es.intime AND es.intime + INTERVAL 24 HOUR
                 THEN rv.val END) AS rr_max_24h,
        MIN(CASE WHEN rv.itemid = 220277 AND rv.val BETWEEN 40.0 AND 100.0
                  AND rv.charttime BETWEEN es.intime AND es.intime + INTERVAL 24 HOUR
                 THEN rv.val END) AS spo2_min_24h,
        MAX(CASE
              WHEN rv.itemid = 223761 AND rv.val BETWEEN 86.0 AND 110.0
                   AND rv.charttime BETWEEN es.intime AND es.intime + INTERVAL 24 HOUR
              THEN (rv.val - 32.0) * 5.0 / 9.0
              WHEN rv.itemid = 223762 AND rv.val BETWEEN 30.0 AND 43.0
                   AND rv.charttime BETWEEN es.intime AND es.intime + INTERVAL 24 HOUR
              THEN rv.val
            END) AS temp_max_24h
    FROM eligible_sepsis es
    LEFT JOIN raw_vitals rv ON es.stay_id = rv.stay_id
    GROUP BY es.stay_id;
    """)

    # Step 6/6: Assemble final analytic dataframe
    print(f"[6/6] Assembling analytic datasets and computing outcomes & indices ({time.time()-t0:.1f}s)...")
    df = con.execute("""
    SELECT
        es.subject_id,
        es.hadm_id,
        es.stay_id,
        es.first_careunit,
        es.age,
        es.gender,
        es.race,
        es.admission_type,
        es.icu_los_days,
        EPOCH(es.dischtime - es.admittime) / 86400.0 AS hosp_los_days,
        es.hospital_expire_flag,
        CASE
            WHEN es.deathtime IS NOT NULL THEN EPOCH(es.deathtime - es.intime) / 86400.0
            WHEN es.dod IS NOT NULL THEN EPOCH((es.dod + INTERVAL 23 HOUR + INTERVAL 59 MINUTE) - es.intime) / 86400.0
            ELSE NULL
        END AS death_days_from_icu,
        es.septic_shock_icd,
        es.chf,
        es.ckd,
        es.liver_disease,
        es.cirrhosis_portal_htn,
        es.diabetes,
        es.malignancy,
        es.copd,
        es.cad_mi,
        es.cerebrovascular,
        es.source_pulmonary,
        es.source_abdominal,
        es.source_urinary,
        es.culture_positive_flag,
        COALESCE(pr.vasopressor_24h, 0) AS vasopressor_24h,
        COALESCE(pr.norepi_24h, 0) AS norepi_24h,
        COALESCE(pr.vasopressin_24h, 0) AS vasopressin_24h,
        COALESCE(pc.imv_24h, 0) AS imv_24h,
        COALESCE(pc.rrt_24h, 0) AS rrt_24h,
        la.lactate_24h,
        la.albumin_24h,
        la.platelet_24h,
        la.inr_24h,
        la.wbc_24h,
        la.creatinine_24h,
        la.bilirubin_24h,
        la.bun_24h,
        la.bicarbonate_24h,
        la.hemoglobin_24h,
        la.lactate_ext24h,
        la.albumin_ext24h,
        la.platelet_ext24h,
        la.inr_ext24h,
        va.map_min_24h,
        va.map_mean_24h,
        va.hr_max_24h,
        va.rr_max_24h,
        va.spo2_min_24h,
        va.temp_max_24h
    FROM eligible_sepsis es
    LEFT JOIN pressors_24h pr ON es.stay_id = pr.stay_id
    LEFT JOIN proc_24h pc ON es.stay_id = pc.stay_id
    LEFT JOIN labs_agg la ON es.stay_id = la.stay_id
    LEFT JOIN vitals_agg va ON es.stay_id = va.stay_id
    ORDER BY es.stay_id
    """).df()

    # Clean demographic categories
    df["sex"] = np.where(df["gender"] == "M", "Male", "Female")
    race_upper = df["race"].fillna("UNKNOWN").str.upper()
    df["race_group"] = np.select(
        [
            race_upper.str.contains("WHITE"),
            race_upper.str.contains("BLACK|AFRICAN"),
            race_upper.str.contains("HISPANIC|LATINO"),
            race_upper.str.contains("ASIAN"),
        ],
        ["White", "Black", "Hispanic", "Asian"],
        default="Other/Unknown",
    )
    adm_upper = df["admission_type"].fillna("UNKNOWN").str.upper()
    df["admission_group"] = np.select(
        [
            adm_upper.str.contains("EMER"),
            adm_upper.str.contains("URGENT|DIRECT"),
        ],
        ["Emergency", "Urgent/Direct"],
        default="Elective/Surgical",
    )

    # Survival endpoints (28-day & 90-day all-cause mortality from ICU admission, post-24h landmark)
    df["mort_28d"] = np.where(df["death_days_from_icu"].notna() & (df["death_days_from_icu"] <= 28.0), 1, 0)
    df["time_28d"] = np.where(
        df["mort_28d"] == 1,
        np.clip(df["death_days_from_icu"], 1.001, 28.0),
        28.0,
    )
    df["time_post24h_27d"] = df["time_28d"] - 1.0

    df["mort_90d"] = np.where(df["death_days_from_icu"].notna() & (df["death_days_from_icu"] <= 90.0), 1, 0)
    df["time_90d"] = np.where(
        df["mort_90d"] == 1,
        np.clip(df["death_days_from_icu"], 1.001, 90.0),
        90.0,
    )

    # Approximate Charlson Comorbidity Burden Score
    age_pts = np.clip((df["age"] - 40) // 10, 0, 4)
    df["charlson_approx"] = (
        age_pts
        + df["chf"]
        + df["cad_mi"]
        + df["cerebrovascular"]
        + df["copd"]
        + df["diabetes"]
        + 2 * df["ckd"]
        + 2 * df["malignancy"]
        + np.where(df["cirrhosis_portal_htn"] == 1, 3, df["liver_disease"])
    )

    # Acute organ dysfunction count (non-exposure clinical severity proxy + SOFA-like burden)
    df["sepsis3_organ_dysfunction"] = (
        (df["vasopressor_24h"] == 1)
        | (df["map_min_24h"].fillna(75) < 70.0)
        | (df["imv_24h"] == 1)
        | (df["creatinine_24h"].fillna(0.8) >= 1.2)
        | (df["bilirubin_24h"].fillna(0.5) >= 1.2)
        | (df["platelet_24h"].fillna(200) < 150.0)
    ).astype(int)

    # Complete-case flag for primary [-6h, +24h] endothelial-coagulopathy triad/tetrad
    df["complete_panel_24h"] = (
        df["lactate_24h"].notna()
        & df["albumin_24h"].notna()
        & df["platelet_24h"].notna()
        & df["inr_24h"].notna()
    ).astype(int)

    # Extended [-24h, +24h] complete-case flag (captures ED pre-ICU admission labs)
    df["complete_panel_ext24h"] = (
        df["lactate_ext24h"].notna()
        & df["albumin_ext24h"].notna()
        & df["platelet_ext24h"].notna()
        & df["inr_ext24h"].notna()
    ).astype(int)

    # Compute LAR and CECI on complete_panel_24h
    # CECI integrates the 4 transcriptomically-anchored axes (standardized z-score composite):
    #   + z(log2(Lactate)) [microvascular hypoperfusion / glycolytic shift]
    #   - z(Albumin)       [endothelial barrier breakdown / capillary leak]
    #   - z(log2(Platelet))[consumptive immunothrombosis]
    #   + z(log2(INR))     [extrinsic/intrinsic coagulopathy]
    df["lar_24h"] = df["lactate_24h"] / df["albumin_24h"]
    df["lar_ext24h"] = df["lactate_ext24h"] / df["albumin_ext24h"]

    # Save full eligible cohort (N = 7,195) for missingness audit & imputation sensitivity
    full_path = PROC_DIR / "mimic_sepsis_eligible_full.csv"
    df.to_csv(full_path, index=False)

    # Primary complete-case cohort
    primary = df[df["complete_panel_24h"] == 1].copy()
    for col, raw_s, sign in [
        ("z_lac", np.log2(primary["lactate_24h"]), +1.0),
        ("z_alb", primary["albumin_24h"], -1.0),
        ("z_plt", np.log2(primary["platelet_24h"]), -1.0),
        ("z_inr", np.log2(primary["inr_24h"]), +1.0),
    ]:
        mu, sd = raw_s.mean(), raw_s.std(ddof=1)
        primary[col] = sign * (raw_s - mu) / sd

    primary["ceci_raw"] = (primary["z_lac"] + primary["z_alb"] + primary["z_plt"] + primary["z_inr"]) / 4.0
    # Standardize CECI to mean=0, SD=1 across primary cohort for direct per-1-SD HR interpretation
    primary["ceci_z"] = (primary["ceci_raw"] - primary["ceci_raw"].mean()) / primary["ceci_raw"].std(ddof=1)

    # Also compute clinical bedside point score (0-8 points, integer clinical translation)
    # Lactate (<2: 0, 2-3.9: 1, >=4: 2), Albumin (>=3.0: 0, 2.5-2.9: 1, <2.5: 2),
    # Platelets (>=150: 0, 100-149: 1, <100: 2), INR (<1.3: 0, 1.3-1.7: 1, >1.7: 2)
    pts_lac = np.select([primary["lactate_24h"] >= 4.0, primary["lactate_24h"] >= 2.0], [2, 1], default=0)
    pts_alb = np.select([primary["albumin_24h"] < 2.5, primary["albumin_24h"] < 3.0], [2, 1], default=0)
    pts_plt = np.select([primary["platelet_24h"] < 100.0, primary["platelet_24h"] < 150.0], [2, 1], default=0)
    pts_inr = np.select([primary["inr_24h"] > 1.7, primary["inr_24h"] >= 1.3], [2, 1], default=0)
    primary["ceci_points"] = pts_lac + pts_alb + pts_plt + pts_inr

    prim_path = PROC_DIR / "mimic_sepsis_endothelial_cohort.csv"
    primary.to_csv(prim_path, index=False)

    # Export STROBE Attrition Table
    strobe_rows = [
        ("Step 1", "All ICU admissions in MIMIC-IV v3.1", int(strobe_pre["step1_all_icu_stays"]), 0),
        ("Step 2", "Adult ICU admissions (age >= 18 years)", int(strobe_pre["step2_adult_stays"]), int(strobe_pre["step1_all_icu_stays"] - strobe_pre["step2_adult_stays"])),
        ("Step 3", "First ICU admission per unique patient", int(strobe_pre["step3_first_icu_stay"]), int(strobe_pre["step2_adult_stays"] - strobe_pre["step3_first_icu_stay"])),
        ("Step 4", "ICU stay >= 24 hours and alive at 24h landmark (zero immortal-time window)", int(strobe_pre["step4_los_ge_24h_landmark"]), int(strobe_pre["step3_first_icu_stay"] - strobe_pre["step4_los_ge_24h_landmark"])),
        ("Step 5", "Confirmed ICU-admission Sepsis (Sepsis ICD-9/10 + active microbiological culture in [-48h, +24h])", int(strobe_pre["step5_sepsis_admission_verified"]), int(strobe_pre["step4_los_ge_24h_landmark"] - strobe_pre["step5_sepsis_admission_verified"])),
        ("Step 6", "Primary Complete-Case Cohort: complete Day-1 [-6h, +24h] Lactate, Albumin, Platelets, and INR", int(len(primary)), int(len(df) - len(primary))),
    ]
    strobe_df = pd.DataFrame(strobe_rows, columns=["step", "criterion", "n_remaining", "n_excluded"])
    strobe_df.to_csv(TAB_DIR / "Table_S1_MIMIC_STROBE_Attrition.csv", index=False)

    # Print summary audit
    print("\n=== MIMIC-IV v3.1 EXTRACTION SUMMARY ===")
    print(strobe_df.to_string(index=False))
    print(f"\nFull Eligible Sepsis Cohort: N = {len(df)}, 28-day deaths = {df['mort_28d'].sum()} ({100*df['mort_28d'].mean():.2f}%)")
    print("Biomarker availability in [-6h, +24h] among eligible sepsis (N=7,195):")
    for col in ["lactate_24h", "albumin_24h", "platelet_24h", "inr_24h", "wbc_24h", "creatinine_24h", "bilirubin_24h", "map_min_24h"]:
        avail = df[col].notna().sum()
        print(f"  {col:16s}: {avail:5d} / {len(df)} ({100*avail/len(df):.1f}%)")
    print(f"Extended [-24h, +24h] complete panel: N = {df['complete_panel_ext24h'].sum()}")
    print(f"Primary [-6h, +24h] Complete-Case Cohort: N = {len(primary)}, 28-day deaths = {primary['mort_28d'].sum()} ({100*primary['mort_28d'].mean():.2f}%)")
    print(f"Completed in {time.time() - t0:.1f}s.")


if __name__ == "__main__":
    main()
