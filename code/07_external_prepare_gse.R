## =========================================================
## 07_external_prepare_gse.R
## Purpose:
##   For chosen external GSE(s): download -> unwrap ExpressionSet ->
##   save pheno/feature info -> generate pheno field report (keywords)
## Outputs (per GSE):
##   01_data_raw/<GSE>/GSE*_eset.rds
##   02_data_clean/<GSE>_pheno.csv
##   02_data_clean/<GSE>_platform.txt
##   00_protocol/<GSE>_pheno_field_report.csv
## =========================================================

source('03_scripts/00_paths.R')
options(stringsAsFactors = FALSE)

if (!requireNamespace('GEOquery', quietly = TRUE)) install.packages('GEOquery')
if (!requireNamespace('Biobase', quietly = TRUE)) install.packages('Biobase')
if (!requireNamespace('dplyr', quietly = TRUE)) install.packages('dplyr')
if (!requireNamespace('stringr', quietly = TRUE)) install.packages('stringr')
if (!requireNamespace('readr', quietly = TRUE)) install.packages('readr')
if (!requireNamespace('tibble', quietly = TRUE)) install.packages('tibble')

suppressPackageStartupMessages({
  library(GEOquery)
  library(Biobase)
  library(dplyr)
  library(stringr)
  library(readr)
  library(tibble)
})

safe_mkdir <- function(d) if (!dir.exists(d)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
timestamp <- function() format(Sys.time(), '%Y-%m-%d %H:%M:%S')

# ---- configure which GSE to prepare ----
# Start with the top-2 you should lock first:
GSE_LIST <- c('GSE26378','GSE28750')
# If you want to also inspect GSE95233, later change to:
# GSE_LIST <- c('GSE26378','GSE28750','GSE95233')

GET_GPL <- TRUE

unwrap_eset <- function(gobj) {
  if (inherits(gobj, 'ExpressionSet')) return(gobj)
  if (is.list(gobj) && length(gobj) >= 1) {
    for (i in seq_along(gobj)) {
      if (inherits(gobj[[i]], 'ExpressionSet')) return(gobj[[i]])
    }
  }
  NULL
}

# Keyword scan across pheno columns, returning per-column hit counts
pheno_field_report <- function(ph) {
  if (is.null(ph) || nrow(ph) == 0) return(tibble())
  keys <- list(
    sepsis = c('sepsis','septic','septic shock'),
    infection = c('infection','infected','suspected'),
    sirs = c('sirs'),
    control = c('control','healthy','volunteer','non-sepsis','nonsepsis'),
    mortality = c('mortality','death','dead','died','surviv','survival'),
    time = c('day','days','time','follow','tte','28')
  )

  cn <- colnames(ph)
  out <- lapply(cn, function(cc) {
    v <- tolower(paste(as.character(ph[[cc]]), collapse=' | '))
    tibble(
      field = cc,
      n_unique = dplyr::n_distinct(ph[[cc]]),
      hit_sepsis = any(stringr::str_detect(v, paste(keys$sepsis, collapse='|'))),
      hit_infection = any(stringr::str_detect(v, paste(keys$infection, collapse='|'))),
      hit_sirs = any(stringr::str_detect(v, paste(keys$sirs, collapse='|'))),
      hit_control = any(stringr::str_detect(v, paste(keys$control, collapse='|'))),
      hit_mortality = any(stringr::str_detect(v, paste(keys$mortality, collapse='|'))),
      hit_time = any(stringr::str_detect(v, paste(keys$time, collapse='|'))),
      example_values = paste(head(unique(as.character(ph[[cc]])), 6), collapse='; ')
    )
  })
  dplyr::bind_rows(out) %>%
    dplyr::arrange(desc(hit_sepsis), desc(hit_mortality), desc(hit_control), desc(hit_infection), desc(hit_time), desc(n_unique))
}

for (gse in GSE_LIST) {
  message('[', timestamp(), '] Downloading ', gse, ' ...')
  gobj <- GEOquery::getGEO(gse, GSEMatrix = TRUE, getGPL = GET_GPL)
  eset <- unwrap_eset(gobj)
  if (is.null(eset)) stop('Cannot unwrap getGEO output to ExpressionSet for: ', gse)

  # Save raw eset for reproducibility
  raw_dir <- file.path(dir_raw, gse)
  safe_mkdir(raw_dir)
  saveRDS(eset, file.path(raw_dir, paste0(gse, '_eset.rds')))

  # Export platform id
  gpl <- tryCatch(Biobase::annotation(eset), error=function(e) NA_character_)
  writeLines(as.character(gpl), con = file.path(dir_clean, paste0(gse, '_platform.txt')))

  # Export pheno
  ph <- Biobase::pData(eset)
  readr::write_csv(as.data.frame(ph), file.path(dir_clean, paste0(gse, '_pheno.csv')))

  # Field report to protocol folder
  rpt <- pheno_field_report(ph)
  readr::write_csv(rpt, file.path(dir_protocol, paste0(gse, '_pheno_field_report.csv')))

  message('[', timestamp(), '] Done: ', gse, ' (n=', ncol(Biobase::exprs(eset)), ', platform=', gpl, ')')
}

cat('✅ External GSE preparation done.\n')
cat('Check outputs in 00_protocol/ and 02_data_clean/.\n')
