# Verification of the code snapshot

- Parsed 10 R and 3 Python scripts.
- Reran the eICU extraction and assessment using an existing local MIMIC-IV development cohort and individually credentialed eICU-CRD v2.0 files. The aggregate assessment CSV matched the manuscript analysis output byte for byte. The complete panel contained 4,169 patients at 159 hospitals, with 984 in-hospital deaths. The fixed CECI AUC was 0.707 (95% hospital-cluster bootstrap CI, 0.688 to 0.727).
- Verified the GEO reference tables against the original validation inputs. The gene order matched exactly and the locked coefficients agreed within 1e-12.
- The full MIMIC-IV extraction and GEO analyses were outside this verification run. The GEO workflow requires the historical preprocessed inputs identified by the input checksums.
- The archive contains source files and documentation. The individually credentialed inputs, clinical cohorts, predictions, and fitted clinical transformations stay in the private work directory.

This note records the checks completed for version 1.0.0.
