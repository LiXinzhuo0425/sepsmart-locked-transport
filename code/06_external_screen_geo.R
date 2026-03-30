# ============================================================
# 06_external_screen_geo.R (v2 FIXED: list-unpack + GPL mapping)
# Purpose:
#   Robust screening of external GEO cohorts for SepSMART validation
# Outputs:
#   04_results/T10_external_candidates_screening.csv
#   04_results/T10_external_candidates_screening_top.tsv
# ============================================================

source('03_scripts/00_paths.R')
options(stringsAsFactors = FALSE)

SEPSMART_GENES <- c('RETN','MCEMP1','CYP1B1','S100A12','S100A8','HK3')

# Starter candidate list (we will expand after the pipeline is correct)
CANDIDATE_GSE <- c('GSE57065','GSE28750','GSE54514','GSE26378','GSE95233','GSE65682')

GET_GPL <- TRUE

bootstrap_pkgs <- function(pkgs) {
  missing <- pkgs[!vapply(pkgs, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))]
  if (length(missing) > 0) {
    message('[PKG] Installing missing packages: ', paste(missing, collapse = ', '))
    if (requireNamespace('renv', quietly = TRUE)) {
      try(renv::install(missing), silent = TRUE)
    }
    missing2 <- missing[!vapply(missing, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))]
    if (length(missing2) > 0) install.packages(missing2)
  }
}

bootstrap_pkgs(c('GEOquery','Biobase','dplyr','tibble','readr','stringr','purrr','data.table'))

suppressPackageStartupMessages({
  library(GEOquery)
  library(Biobase)
  library(dplyr)
  library(tibble)
  library(readr)
  library(stringr)
  library(purrr)
  library(data.table)
})

timestamp <- function() format(Sys.time(), '%Y-%m-%d %H:%M:%S')

norm_symbol_vec <- function(x) {
  x <- as.character(x); x[is.na(x)] <- ''
  x <- toupper(x)
  x <- str_replace_all(x, '///', ';')
  x <- str_replace_all(x, '\\|', ';')
  x <- str_replace_all(x, ',', ';')
  x <- str_replace_all(x, '\\s+', '')
  x
}

# Robustly unwrap getGEO output to an ExpressionSet
unwrap_eset <- function(gobj) {
  if (inherits(gobj, 'ExpressionSet')) return(gobj)
  if (is.list(gobj) && length(gobj) >= 1) {
    # many GSE return list(ExpressionSet) even when length==1
    for (i in seq_along(gobj)) {
      if (inherits(gobj[[i]], 'ExpressionSet')) return(gobj[[i]])
    }
  }
  NULL
}

detect_pheno_signals <- function(pheno_df) {
  if (is.null(pheno_df) || nrow(pheno_df) == 0) {
    return(list(has_sepsis=NA, has_sirs=NA, has_infection=NA, has_control=NA, has_mortality=NA, has_days=NA))
  }
  txt <- apply(pheno_df, 2, function(col) paste(as.character(col), collapse = ' | '))
  t <- tolower(paste(txt, collapse = ' || '))
  list(
    has_sepsis = str_detect(t, '\\bsepsis\\b|septic shock|septic'),
    has_sirs = str_detect(t, '\\bsirs\\b'),
    has_infection = str_detect(t, 'infection|infected|suspected'),
    has_control = str_detect(t, 'control|healthy|volunteer|non-sepsis|nonsepsis'),
    has_mortality = str_detect(t, 'mortality|death|dead|surviv|survival|died'),
    has_days = str_detect(t, 'day|days|time-to|time_to|follow|28')
  )
}

# Try to extract gene symbols from fData first
extract_gene_symbols_from_fdata <- function(eset) {
  fd <- tryCatch(Biobase::fData(eset), error = function(e) NULL)
  if (is.null(fd) || nrow(fd) == 0) return(character(0))
  cn <- colnames(fd)
  sym_cols <- cn[str_detect(tolower(cn), 'gene.?symbol|symbol|genesymbol')]
  if (length(sym_cols) == 0) return(character(0))
  scores <- sapply(sym_cols, function(cc) sum(!is.na(fd[[cc]]) & fd[[cc]] != ''))
  sym_col <- sym_cols[which.max(scores)]
  syms <- norm_symbol_vec(fd[[sym_col]])
  syms <- unlist(str_split(syms, ';'), use.names = FALSE)
  syms <- syms[syms != '']
  unique(syms)
}

# If fData lacks gene symbols, use GPL table to map probes -> symbols
probe_to_symbol_via_gpl <- function(eset, gpl_id) {
  probes <- rownames(Biobase::exprs(eset))
  if (is.null(probes) || length(probes) == 0) return(list(symbols=character(0), note='No probes'))

  gpl_obj <- tryCatch(GEOquery::getGEO(gpl_id), error = function(e) NULL)
  if (is.null(gpl_obj)) return(list(symbols=character(0), note='Failed to download GPL'))

  tab <- tryCatch(GEOquery::Table(gpl_obj), error = function(e) NULL)
  if (is.null(tab) || nrow(tab) == 0) return(list(symbols=character(0), note='GPL table empty'))

  # Identify probe id column
  id_col <- if ('ID' %in% colnames(tab)) 'ID' else colnames(tab)[1]

  # Identify symbol column candidates
  cn <- colnames(tab)
  sym_cols <- cn[str_detect(tolower(cn), 'gene.?symbol|symbol|gene symbol')]
  if (length(sym_cols) == 0) {
    return(list(symbols=character(0), note='No symbol-like column in GPL table'))
  }
  scores <- sapply(sym_cols, function(cc) sum(!is.na(tab[[cc]]) & tab[[cc]] != ''))
  sym_col <- sym_cols[which.max(scores)]

  tab2 <- tab[, c(id_col, sym_col)]
  colnames(tab2) <- c('probe','symbol')
  tab2$symbol <- norm_symbol_vec(tab2$symbol)
  tab2 <- tab2[tab2$probe %in% probes, , drop=FALSE]
  if (nrow(tab2) == 0) return(list(symbols=character(0), note='No probe overlap between eset and GPL table'))

  syms <- unlist(str_split(tab2$symbol, ';'), use.names = FALSE)
  syms <- syms[syms != '']
  list(symbols=unique(syms), note=paste0('Mapped via ', gpl_id, ' / ', sym_col))
}

safe_getGEO <- function(gse_id, get_gpl = TRUE) {
  message('[', timestamp(), '] Fetching ', gse_id, ' ...')
  out <- tryCatch({
    gset <- GEOquery::getGEO(gse_id, GSEMatrix = TRUE, getGPL = get_gpl)
    gset
  }, error = function(e) {
    structure(list(error = TRUE, message = as.character(e$message)), class = 'geo_error')
  })
  out
}

results <- purrr::map_dfr(CANDIDATE_GSE, function(gse) {
  gobj <- safe_getGEO(gse, get_gpl = GET_GPL)
  if (inherits(gobj, 'geo_error')) {
    return(tibble(gse=gse, status='FAIL_DOWNLOAD', n_samples=NA_integer_, platform=NA_character_,
      has_sepsis=NA, has_sirs=NA, has_infection=NA, has_control=NA, has_mortality=NA, has_days=NA,
      genes_present_n=NA_integer_, genes_missing=NA_character_, note=gobj$message))
  }

  eset <- unwrap_eset(gobj)
  if (is.null(eset)) {
    return(tibble(gse=gse, status='FAIL_PARSE', n_samples=NA_integer_, platform=NA_character_,
      has_sepsis=NA, has_sirs=NA, has_infection=NA, has_control=NA, has_mortality=NA, has_days=NA,
      genes_present_n=NA_integer_, genes_missing=paste(SEPSMART_GENES, collapse=';'),
      note='Cannot unwrap to ExpressionSet (unexpected getGEO object)'))
  }

  n_samples <- ncol(Biobase::exprs(eset))
  gpl_id <- tryCatch(Biobase::annotation(eset), error=function(e) NA_character_)

  ph <- tryCatch(Biobase::pData(eset), error=function(e) data.frame())
  sig <- detect_pheno_signals(ph)

  # Gene presence: fData symbols first, else GPL mapping
  syms <- extract_gene_symbols_from_fdata(eset)
  note_gene <- ''
  if (length(syms) == 0 && !is.na(gpl_id) && nzchar(gpl_id)) {
    mapped <- probe_to_symbol_via_gpl(eset, gpl_id)
    syms <- mapped$symbols
    note_gene <- mapped$note
  } else if (length(syms) > 0) {
    note_gene <- 'Symbols from fData'
  }

  if (length(syms) == 0) {
    present_n <- NA_integer_
    missing <- paste(SEPSMART_GENES, collapse=';')
    note_gene2 <- ifelse(nzchar(note_gene), note_gene, 'No gene symbols after mapping')
  } else {
    present <- SEPSMART_GENES %in% syms
    present_n <- sum(present)
    missing <- paste(SEPSMART_GENES[!present], collapse=';')
    note_gene2 <- paste0(note_gene, '; present=', present_n, '/6')
  }

  tibble(
    gse=gse, status='OK', n_samples=n_samples, platform=gpl_id,
    has_sepsis=sig$has_sepsis, has_sirs=sig$has_sirs, has_infection=sig$has_infection, has_control=sig$has_control,
    has_mortality=sig$has_mortality, has_days=sig$has_days,
    genes_present_n=present_n, genes_missing=missing, note=note_gene2
  )
})

ranked <- results %>%
  mutate(
    ok = status=='OK',
    genes_present_n2 = ifelse(is.na(genes_present_n), 0, genes_present_n),
    score = 100*ok + 10*genes_present_n2 +
      5*(has_sepsis %in% TRUE) + 3*(has_control %in% TRUE) + 2*(has_infection %in% TRUE) +
      2*(has_sirs %in% TRUE) + 2*(has_mortality %in% TRUE) + 1*(has_days %in% TRUE)
  ) %>%
  arrange(desc(score), desc(genes_present_n2), desc(n_samples))

out_csv <- file.path(dir_results, 'T10_external_candidates_screening.csv')
out_top <- file.path(dir_results, 'T10_external_candidates_screening_top.tsv')
readr::write_csv(ranked, out_csv)
readr::write_tsv(ranked %>% slice(1:min(20, n())), out_top)

cat('✅ External candidate screening done.\n')
cat(' - Full table: ', out_csv, '\n', sep='')
cat(' - Top preview: ', out_top, '\n', sep='')
