#!/usr/bin/env python3
"""
05_extract_dynamic_and_fluid_tte.py
===================================
Refactored and mathematically harmonized extraction of:
1. True 24-hour crystalloid resuscitation volume from inputevents.csv:
   - Window strictly constrained to [intime, intime + 24h]
   - Unit conversion (L -> mL)
   - Time-proportional allocation for cross-boundary infusions
   - Excludes rewritten events
2. Serial Day 2 (24-48h) and Day 3 (48-72h) labs from labevents.csv
3. Longitudinal CECI computation strictly locked to baseline formulas & standardization parameters:
   - log2(lactate), albumin, log2(platelets), log2(INR)
   - Baseline mean & standard deviation parameters locked from Day 1
   - True within-person delta = ceci_d2_z - ceci_z
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

def main():
    t0 = time.time()
    print("=== [1/4] Loading Primary MIMIC Sepsis Cohort ===")
    cohort = pd.read_csv(PROC_DIR / "mimic_sepsis_endothelial_cohort.csv")
    print(f"Loaded primary cohort: {len(cohort)} ICU stays")
    
    con = duckdb.connect(":memory:")
    con.execute("SET threads=6; SET memory_limit='6GB'")
    con.register("cohort_df", cohort[["stay_id", "hadm_id", "subject_id"]])
    
    print("\n=== [2/4] Extracting True Windowed 24h Crystalloid Volume from inputevents.csv ===")
    fluid_q = f"""
    WITH stays AS (
        SELECT c.stay_id, i.intime, i.intime + INTERVAL 24 HOUR AS end_24h
        FROM cohort_df c
        JOIN read_csv_auto('{MIMIC_DIR.as_posix()}/icu/icustays.csv') i ON c.stay_id = i.stay_id
    ),
    events AS (
        SELECT 
            s.stay_id,
            TRY_CAST(ie.patientweight AS DOUBLE) as patientweight,
            TRY_CAST(ie.starttime AS TIMESTAMP) as st,
            COALESCE(TRY_CAST(ie.endtime AS TIMESTAMP), TRY_CAST(ie.starttime AS TIMESTAMP)) as et,
            CASE WHEN LOWER(TRIM(ie.amountuom)) = 'l' THEN TRY_CAST(ie.amount AS DOUBLE) * 1000.0 ELSE TRY_CAST(ie.amount AS DOUBLE) END as amt_ml,
            ie.itemid,
            ie.statusdescription,
            s.intime,
            s.end_24h
        FROM stays s
        JOIN read_csv_auto('{MIMIC_DIR.as_posix()}/icu/inputevents.csv') ie ON s.stay_id = ie.stay_id
        WHERE ie.itemid IN (225158, 225828, 225823, 225825, 225827, 225159, 226364, 226375, 220949, 225161, 220862, 220864, 220861, 220954, 220955, 220956, 221208)
          AND COALESCE(ie.statusdescription, '') != 'Rewritten'
          AND LOWER(TRIM(ie.amountuom)) IN ('ml', 'l')
    ),
    windowed AS (
        SELECT 
            stay_id,
            itemid,
            patientweight,
            CASE 
                WHEN amt_ml IS NULL OR amt_ml <= 0 THEN 0.0
                WHEN st >= end_24h OR et < intime THEN 0.0
                WHEN et <= st THEN 
                    CASE WHEN st >= intime AND st < end_24h THEN amt_ml ELSE 0.0 END
                ELSE 
                    amt_ml * (EPOCH(LEAST(et, end_24h)) - EPOCH(GREATEST(st, intime))) / (EPOCH(et) - EPOCH(st))
            END AS vol_ml_24h
        FROM events
    )
    SELECT 
        c.stay_id,
        COALESCE(AVG(CASE WHEN w.vol_ml_24h > 0 AND w.patientweight BETWEEN 40 AND 180 THEN w.patientweight END), 75.0) as weight_kg,
        COALESCE(SUM(CASE WHEN w.itemid IN (225828, 220955, 220956, 221208, 225158, 220954) THEN w.vol_ml_24h ELSE 0 END), 0.0) as crystalloid_vol_24h,
        COALESCE(SUM(CASE WHEN w.itemid IN (225828, 220955, 220956, 221208) THEN w.vol_ml_24h ELSE 0 END), 0.0) as balanced_vol_24h,
        COALESCE(SUM(CASE WHEN w.itemid = 225827 THEN w.vol_ml_24h ELSE 0 END), 0.0) as d5lr_vol_24h,
        COALESCE(SUM(CASE WHEN w.itemid IN (225158, 220954) THEN w.vol_ml_24h ELSE 0 END), 0.0) as saline_vol_24h,
        COALESCE(SUM(CASE WHEN w.itemid IN (220862, 220864, 220861) THEN w.vol_ml_24h ELSE 0 END), 0.0) as albumin_vol_24h
    FROM cohort_df c
    LEFT JOIN windowed w ON c.stay_id = w.stay_id
    GROUP BY c.stay_id
    """
    fluid_df = con.execute(fluid_q).df()
    print(f"Extracted fluid data for {len(fluid_df)} stays.")
    print("Fluid volume summary (mL):")
    print(fluid_df["crystalloid_vol_24h"].describe())
    
    print("\n=== [3/4] Extracting TRUE Day-2 (24-48h) & Day-3 (48-72h) Biomarkers ===")
    lab_q = f"""
    CREATE TEMP TABLE cand_hadm AS SELECT DISTINCT hadm_id FROM cohort_df;
    
    CREATE TEMP TABLE raw_serial_labs AS
    SELECT 
        l.hadm_id,
        TRY_CAST(l.charttime AS TIMESTAMP) AS charttime,
        l.itemid,
        TRY_CAST(l.valuenum AS DOUBLE) AS valuenum
    FROM read_csv('{MIMIC_DIR.as_posix()}/hosp/labevents.csv', header=true, all_varchar=true, sample_size=20480) l
    JOIN cand_hadm c ON l.hadm_id = c.hadm_id
    WHERE l.itemid IN ('50813', '52442', '53154', '50862', '53085', '53138', '52022', '51265', '53189', '51237', '51675', '50912')
      AND l.valuenum IS NOT NULL;

    CREATE TEMP TABLE stay_windows AS
    SELECT 
        c.stay_id,
        c.hadm_id,
        i.intime,
        i.intime + INTERVAL 24 HOUR AS t24,
        i.intime + INTERVAL 48 HOUR AS t48,
        i.intime + INTERVAL 72 HOUR AS t72, h.dischtime
    FROM cohort_df c
    JOIN read_csv_auto('{MIMIC_DIR.as_posix()}/icu/icustays.csv') i ON c.stay_id = i.stay_id
    JOIN read_csv_auto('{MIMIC_DIR.as_posix()}/hosp/admissions.csv') h ON h.hadm_id = c.hadm_id;

    SELECT 
        sw.stay_id,
        -- Day 2: 24h to 48h (Worst physiological value)
        MAX(CASE WHEN r.itemid IN ('50813', '52442', '53154') AND r.valuenum BETWEEN 0.1 AND 30 AND r.charttime >= sw.t24 AND r.charttime < sw.t48 THEN r.valuenum END) AS lactate_d2,
        MIN(CASE WHEN r.itemid IN ('50862', '53085', '53138', '52022') AND r.valuenum BETWEEN 0.5 AND 6.5 AND r.charttime >= sw.t24 AND r.charttime < sw.t48 THEN r.valuenum END) AS albumin_d2,
        MIN(CASE WHEN r.itemid IN ('51265', '53189') AND r.valuenum BETWEEN 1 AND 1500 AND r.charttime >= sw.t24 AND r.charttime < sw.t48 THEN r.valuenum END) AS platelet_d2,
        MAX(CASE WHEN r.itemid IN ('51237', '51675') AND r.valuenum BETWEEN 0.5 AND 25 AND r.charttime >= sw.t24 AND r.charttime < sw.t48 THEN r.valuenum END) AS inr_d2,
        -- Day 3: 48h to 72h
        MAX(CASE WHEN r.itemid IN ('50813', '52442', '53154') AND r.valuenum BETWEEN 0.1 AND 30 AND r.charttime >= sw.t48 AND r.charttime < sw.t72 THEN r.valuenum END) AS lactate_d3,
        MIN(CASE WHEN r.itemid IN ('50862', '53085', '53138', '52022') AND r.valuenum BETWEEN 0.5 AND 6.5 AND r.charttime >= sw.t48 AND r.charttime < sw.t72 THEN r.valuenum END) AS albumin_d3,
        MIN(CASE WHEN r.itemid IN ('51265', '53189') AND r.valuenum BETWEEN 1 AND 1500 AND r.charttime >= sw.t48 AND r.charttime < sw.t72 THEN r.valuenum END) AS platelet_d3,
        MAX(CASE WHEN r.itemid IN ('51237', '51675') AND r.valuenum BETWEEN 0.5 AND 25 AND r.charttime >= sw.t48 AND r.charttime < sw.t72 THEN r.valuenum END) AS inr_d3
    FROM stay_windows sw
    LEFT JOIN raw_serial_labs r ON sw.hadm_id = r.hadm_id AND r.charttime <= sw.dischtime
    GROUP BY sw.stay_id
    """
    serial_labs_df = con.execute(lab_q).df()
    print(f"Extracted serial labs for {len(serial_labs_df)} stays.")
    print("Non-null D2 counts:", serial_labs_df[["lactate_d2", "albumin_d2", "platelet_d2", "inr_d2"]].notna().sum().to_dict())
    print("Non-null D3 counts:", serial_labs_df[["lactate_d3", "albumin_d3", "platelet_d3", "inr_d3"]].notna().sum().to_dict())
    
    print("\n=== [4/4] Merging and Constructing Harmonized Trajectories & TTE Dataset ===")
    merged = cohort.merge(fluid_df, on="stay_id", how="left")
    merged = merged.merge(serial_labs_df, on="stay_id", how="left")
    
    # Fluid dose calculation (mL/kg)
    merged["weight_kg_clamped"] = merged["weight_kg"].clip(lower=40.0, upper=180.0)
    merged["fluid_ml_per_kg_24h"] = merged["crystalloid_vol_24h"] / merged["weight_kg_clamped"]
    
    # Target Trial Emulation (TTE) Treatment Assignment:
    merged["tte_fluid_group"] = np.where(merged["crystalloid_vol_24h"] > 0, np.where(merged["fluid_ml_per_kg_24h"] >= 30.0, "Liberal_ge_30ml_kg", "Restrictive_lt_30ml_kg"), None)
    merged["tte_fluid_treat"] = np.where(merged["crystalloid_vol_24h"] > 0, (merged["fluid_ml_per_kg_24h"] >= 30.0).astype(int), np.nan)
    
    # Fluid composition & Albumin definitions
    total_cryst = merged["balanced_vol_24h"] + merged["saline_vol_24h"]
    merged["pct_balanced"] = np.where(total_cryst > 0, merged["balanced_vol_24h"] / total_cryst, np.nan)
    merged["predom_balanced"] = np.where(total_cryst > 0, (merged["pct_balanced"] >= 0.5).astype(int), np.nan)
    merged["early_albumin"] = np.where(merged["albumin_vol_24h"] > 0, 1, 0)
    
    # Harmonized Standardization: EXACT BASELINE FORMULA & PARAMETERS
    # log2(lactate) (+), albumin (-), log2(platelets) (-), log2(INR) (+)
    mu_lac = np.log2(cohort["lactate_24h"]).mean()
    sd_lac = np.log2(cohort["lactate_24h"]).std(ddof=1)
    
    mu_alb = cohort["albumin_24h"].mean()
    sd_alb = cohort["albumin_24h"].std(ddof=1)
    
    mu_plt = np.log2(cohort["platelet_24h"]).mean()
    sd_plt = np.log2(cohort["platelet_24h"]).std(ddof=1)
    
    mu_inr = np.log2(cohort["inr_24h"]).mean()
    sd_inr = np.log2(cohort["inr_24h"]).std(ddof=1)
    
    # Verify baseline ceci_raw mean and std
    b_lac = +1.0 * (np.log2(cohort["lactate_24h"]) - mu_lac) / sd_lac
    b_alb = -1.0 * (cohort["albumin_24h"] - mu_alb) / sd_alb
    b_plt = -1.0 * (np.log2(cohort["platelet_24h"]) - mu_plt) / sd_plt
    b_inr = +1.0 * (np.log2(cohort["inr_24h"]) - mu_inr) / sd_inr
    b_raw = (b_lac + b_alb + b_plt + b_inr) / 4.0
    mu_ceci_raw = b_raw.mean()
    sd_ceci_raw = b_raw.std(ddof=1)
    
    # Impute missing D2 with D1 if test unrepeated or discharged
    lac_d2 = merged["lactate_d2"].fillna(merged["lactate_24h"])
    alb_d2 = merged["albumin_d2"].fillna(merged["albumin_24h"])
    plt_d2 = merged["platelet_d2"].fillna(merged["platelet_24h"])
    inr_d2 = merged["inr_d2"].fillna(merged["inr_24h"])
    
    z_lac_d2 = +1.0 * (np.log2(lac_d2) - mu_lac) / sd_lac
    z_alb_d2 = -1.0 * (alb_d2 - mu_alb) / sd_alb
    z_plt_d2 = -1.0 * (np.log2(plt_d2) - mu_plt) / sd_plt
    z_inr_d2 = +1.0 * (np.log2(inr_d2) - mu_inr) / sd_inr
    ceci_d2_raw = (z_lac_d2 + z_alb_d2 + z_plt_d2 + z_inr_d2) / 4.0
    merged["ceci_d2_z"] = (ceci_d2_raw - mu_ceci_raw) / sd_ceci_raw
    
    # Day 3
    lac_d3 = merged["lactate_d3"].fillna(lac_d2)
    alb_d3 = merged["albumin_d3"].fillna(alb_d2)
    plt_d3 = merged["platelet_d3"].fillna(plt_d2)
    inr_d3 = merged["inr_d3"].fillna(inr_d2)
    
    z_lac_d3 = +1.0 * (np.log2(lac_d3) - mu_lac) / sd_lac
    z_alb_d3 = -1.0 * (alb_d3 - mu_alb) / sd_alb
    z_plt_d3 = -1.0 * (np.log2(plt_d3) - mu_plt) / sd_plt
    z_inr_d3 = +1.0 * (np.log2(inr_d3) - mu_inr) / sd_inr
    ceci_d3_raw = (z_lac_d3 + z_alb_d3 + z_plt_d3 + z_inr_d3) / 4.0
    merged["ceci_d3_z"] = (ceci_d3_raw - mu_ceci_raw) / sd_ceci_raw
    
    # True Within-Person Delta CECI (Day 2 - Day 1) on identical standard deviation scale
    merged["delta_ceci_d2_d1"] = merged["ceci_d2_z"] - merged["ceci_z"]
    
    # Harmonized Trajectory Classification:
    # delta <= -0.25 SD: Resolving Injury
    # delta >= +0.25 SD: Deteriorating Injury
    # delta between -0.25 and +0.25:
    #   if baseline ceci_z >= 0: Persistent High Injury
    #   else: Persistent Low Injury
    def classify_traj(row):
        d = row["delta_ceci_d2_d1"]
        b = row["ceci_z"]
        if d <= -0.25:
            return "Resolving_Injury"
        elif d >= 0.25:
            return "Deteriorating_Injury"
        elif b >= 0.0:
            return "Persistent_High_Injury"
        else:
            return "Persistent_Low_Injury"
            
    merged["ceci_trajectory_group"] = merged.apply(classify_traj, axis=1)
    
    out_file = PROC_DIR / "mimic_sepsis_tte_trajectory_cohort.csv"
    merged.to_csv(out_file, index=False)
    print(f"\n[SUCCESS] Saved harmonized serial TTE & Trajectory cohort: {out_file} ({len(merged)} rows)")
    print("\nHarmonized trajectory distribution across full cohort (N=3,254):")
    print(merged["ceci_trajectory_group"].value_counts())
    
    surv_48h = merged[merged["time_28d"] > 2.0]
    print(f"\nHarmonized trajectory distribution in 48h landmark survivors (N={len(surv_48h)}):")
    print(surv_48h["ceci_trajectory_group"].value_counts())
    
    print("\nFluid resuscitation distribution:")
    print(merged["tte_fluid_group"].value_counts())
    print("\nFluid volume summary by treatment group:")
    print(merged.groupby("tte_fluid_group")[["crystalloid_vol_24h", "fluid_ml_per_kg_24h"]].agg(["mean", "std", "median"]))
    print(f"\nTotal elapsed: {time.time() - t0:.1f}s")

if __name__ == "__main__":
    main()
