## =========================================================
## 03_label_map_GSE65682.R
## Build analysis-ready subsets from GSE65682 phenotype table
## Outputs:
##  - 02_data_clean/GSE65682_abdo_casecontrol_{expr,meta}.rds (+csv)
##  - 02_data_clean/GSE65682_mortality28d_{expr,meta}.rds (+csv)
##  - 04_results/T0_GSE65682_flow_counts.csv
## =========================================================

source(file.path('03_scripts', '00_paths.R'))
options(renv.consent = TRUE)

suppressPackageStartupMessages({
  library(data.table)
})

safe_mkdir <- function(d) if (!dir.exists(d)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
save_rds <- function(x, path) { safe_mkdir(dirname(path)); saveRDS(x, path); msg('Saved: ', path) }
write_csv <- function(df, path) { safe_mkdir(dirname(path)); data.table::fwrite(df, path); msg('Saved: ', path) }

## ---- locate inputs ----
expr_path <- list.files(dir_clean, pattern = '^GSE65682_.*_expr_clean\\.rds$', full.names = TRUE)
pheno_path <- list.files(dir_clean, pattern = '^GSE65682_.*_pheno\\.csv$', full.names = TRUE)

if (length(expr_path) == 0) stop('Cannot find expr_clean.rds for GSE65682 under: ', dir_clean)
if (length(pheno_path) == 0) stop('Cannot find pheno.csv for GSE65682 under: ', dir_clean)

# If multiple matches, take the first (should be single_platform)
expr_path <- expr_path[1]
pheno_path <- pheno_path[1]

msg('Loading expr: ', expr_path)
expr <- readRDS(expr_path)
msg('Loading pheno: ', pheno_path)
pheno <- data.table::fread(pheno_path, data.table = FALSE)

## Standardize sample IDs
# In our pipeline, columns(expr) are sample IDs; in pheno CSV, usually rownames were lost to a column named 'rownames' or similar.
candidate_id_cols <- c('geo_accession','Geo_accession','geo_accession_id','accession','sample','Sample','title')
id_col <- intersect(candidate_id_cols, colnames(pheno))
if (length(id_col) == 0) {
  # fallback: try to find any column containing 'GSM'
  gsm_cols <- which(sapply(pheno, function(x) any(grepl('^GSM', as.character(x)))))
  if (length(gsm_cols) > 0) id_col <- colnames(pheno)[gsm_cols[1]]
}
if (length(id_col) == 0) stop('Cannot identify sample ID column in pheno. Please inspect pheno columns.')
id_col <- id_col[1]
pheno$.sample_id <- as.character(pheno[[id_col]])

## Align samples between expr and pheno
expr_samples <- colnames(expr)
pheno_samples <- pheno$.sample_id
common <- intersect(expr_samples, pheno_samples)

flow <- data.frame(step = character(), n = integer(), stringsAsFactors = FALSE)
flow <- rbind(flow, data.frame(step = 'Expr samples (raw)', n = length(expr_samples)))
flow <- rbind(flow, data.frame(step = 'Pheno rows (raw)', n = nrow(pheno)))
flow <- rbind(flow, data.frame(step = 'Common samples (expr ∩ pheno)', n = length(common)))

if (length(common) < 50) stop('Too few common samples after alignment: ', length(common))

expr_aligned <- expr[, common, drop = FALSE]
pheno_aligned <- pheno[match(common, pheno$.sample_id), , drop = FALSE]

## ---- subset A: abdominal sepsis vs GI control ----
label_col_A <- 'abdominal_sepsis_and_controls:ch1'
if (!label_col_A %in% colnames(pheno_aligned)) {
  stop('Missing label column for subset A: ', label_col_A)
}

A_raw <- as.character(pheno_aligned[[label_col_A]])
keep_A <- !is.na(A_raw) & A_raw %in% c('abdo_s','ctrl_GI')
flow <- rbind(flow, data.frame(step = 'Subset A eligible (abdo_s/ctrl_GI non-NA)', n = sum(keep_A)))

expr_A <- expr_aligned[, keep_A, drop = FALSE]
meta_A <- pheno_aligned[keep_A, , drop = FALSE]

meta_A$y_case <- ifelse(as.character(meta_A[[label_col_A]]) == 'abdo_s', 1L, 0L)
meta_A$endpoint <- 'abdominal_sepsis_casecontrol'
meta_A$label_source <- label_col_A

flow <- rbind(flow, data.frame(step = 'Subset A cases (abdo_s)', n = sum(meta_A$y_case == 1L)))
flow <- rbind(flow, data.frame(step = 'Subset A controls (ctrl_GI)', n = sum(meta_A$y_case == 0L)))

## Save subset A
save_rds(expr_A, file.path(dir_clean, 'GSE65682_abdo_casecontrol_expr.rds'))
save_rds(meta_A, file.path(dir_clean, 'GSE65682_abdo_casecontrol_meta.rds'))
write_csv(data.frame(sample = colnames(expr_A), y_case = meta_A$y_case),
          file.path(dir_clean, 'GSE65682_abdo_casecontrol_samples.csv'))

## ---- subset B: 28-day mortality ----
label_col_B <- 'mortality_event_28days:ch1'
time_col_B  <- 'time_to_event_28days:ch1'

if (!label_col_B %in% colnames(pheno_aligned)) stop('Missing label column for subset B: ', label_col_B)

B_raw <- as.character(pheno_aligned[[label_col_B]])
keep_B <- !is.na(B_raw) & B_raw %in% c('0.0','1.0', '0', '1', 0, 1)
flow <- rbind(flow, data.frame(step = 'Subset B eligible (mortality 0/1 non-NA)', n = sum(keep_B)))

expr_B <- expr_aligned[, keep_B, drop = FALSE]
meta_B <- pheno_aligned[keep_B, , drop = FALSE]

meta_B$y_death28 <- as.integer(as.character(meta_B[[label_col_B]]))
meta_B$endpoint <- 'mortality_28days'
meta_B$label_source <- label_col_B

flow <- rbind(flow, data.frame(step = 'Subset B deaths (1)', n = sum(meta_B$y_death28 == 1L, na.rm = TRUE)))
flow <- rbind(flow, data.frame(step = 'Subset B survivors (0)', n = sum(meta_B$y_death28 == 0L, na.rm = TRUE)))

# Optional: attach time-to-event if available
if (time_col_B %in% colnames(meta_B)) {
  meta_B$tte28 <- suppressWarnings(as.numeric(as.character(meta_B[[time_col_B]])))
  msg('Subset B: time_to_event column detected; non-NA=', sum(!is.na(meta_B$tte28)))
} else {
  meta_B$tte28 <- NA_real_
  msg('Subset B: time_to_event column not found; tte28 set to NA.')
}

## Save subset B
save_rds(expr_B, file.path(dir_clean, 'GSE65682_mortality28d_expr.rds'))
save_rds(meta_B, file.path(dir_clean, 'GSE65682_mortality28d_meta.rds'))
write_csv(data.frame(sample = colnames(expr_B), y_death28 = meta_B$y_death28, tte28 = meta_B$tte28),
          file.path(dir_clean, 'GSE65682_mortality28d_samples.csv'))

## ---- save flow table ----
write_csv(flow, file.path(dir_results, 'T0_GSE65682_flow_counts.csv'))

## ---- report SepSMART gene availability in aligned expression ----
have <- intersect(genes_sepSMART, rownames(expr_aligned))
miss <- setdiff(genes_sepSMART, rownames(expr_aligned))
gene_presence <- data.frame(gene = genes_sepSMART, present = genes_sepSMART %in% rownames(expr_aligned))
write_csv(gene_presence, file.path(dir_results, 'T0_GSE65682_SepSMART_gene_presence_aligned.csv'))
msg('SepSMART present: ', paste(have, collapse = ', '))
if (length(miss) > 0) msg('SepSMART missing: ', paste(miss, collapse = ', '))

msg('Done. Created two analysis-ready subsets under 02_data_clean/ and flow table under 04_results/.')
