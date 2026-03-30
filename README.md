# sepsmart-locked-transport

Code repository for archived GEO-based locked-model development and external transport evaluation of a six-gene sepsis host-response signature.

## Overview

This repository contains the analysis code used for archived cohort construction, platform-specific probe-to-gene mapping, six-gene input generation, locked-model development, external transport evaluation, threshold-related analyses, calibration analyses, and generation of selected supplementary outputs.

The repository is intended to reproduce the archived GEO-based workflow reported in the manuscript. It does **not** implement raw CEL-level reprocessing or empirical cross-platform batch correction.

## Study scope

The analytical workflow in this repository was designed to evaluate transport behavior of a locked six-gene host-response model across archived heterogeneous GEO cohorts, rather than to create a newly harmonized multi-cohort expression compendium.

## Core retained GEO accessions

- GSE65682
- GSE28750
- GSE26378

Original archived records should be retrieved directly from GEO using the accession identifiers reported in the manuscript and supplementary materials.

## Repository structure

- `code/`  
  Core analysis scripts for environment checking, GEO download/cleaning, probe-to-gene mapping, locked-model training, external cohort preparation, transport evaluation, threshold analyses, calibration analyses, and audit/provenance checks.

- `results_manifest/`  
  Small manifest-style outputs supporting reproducibility of supplementary tables and provenance reporting.

## Main scripts

- `01_check_env.R`
- `02_geo_download_clean.R`
- `03_label_map_GSE65682.R`
- `06_external_screen_geo.R`
- `07_external_prepare_gse.R`
- `08_external_define_groups.R`
- `09_external_extract_and_roc_REDO.R`
- `10_external_metrics_and_thresholds.R`
- `11_train_fixed_dx_model_on_GSE65682.R`
- `13_fixed_model_threshold_metrics.R`
- `14_external_calibration_and_thresholding.R`
- `99_audit_cv_and_probe_aggregation.R`

## Reproducibility notes

- Archived GEO expression objects were used as downloaded.
- No raw CEL-level reprocessing was performed in this repository workflow.
- No empirical cross-platform renormalization or batch correction was applied.
- Platform-specific probe-to-gene mapping and analysis-stage construction of six-gene inputs are documented in the manuscript supplementary materials.

## Supplementary provenance linkage

- Supplementary Table S6 documents final platform-specific probe manifests and construction rules for six-gene input generation.
- Supplementary Table S7 documents GEO search strategy and accession-level screening.
- Supplementary Table S8 documents archived expression-object provenance, preprocessing metadata availability, and analysis-stage handling.

## Requirements

This repository uses R-based analysis scripts. Package requirements should be reconstructed from the individual scripts and local session records. A more explicit dependency record may be added in a later release.

## Data availability note

This repository does not redistribute large archived GEO expression objects or third-party datasets beyond small manifest-style files. Users should retrieve the original archived records directly from GEO.

## License

License information will be added before the citable release is finalized.

