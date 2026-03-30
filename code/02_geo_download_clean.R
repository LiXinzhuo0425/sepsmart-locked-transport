## =========================================================
## 02_geo_download_clean.R
## Download GEO dataset -> clean expression + pheno -> QC + 6-gene check
## Usage (in console):
##   source('03_scripts/02_geo_download_clean.R')
##   run_geo_pipeline('GSE65682', gpl_pick = NULL)
## =========================================================

source(file.path('03_scripts', '00_paths.R'))
options(renv.consent = TRUE)

suppressPackageStartupMessages({
  library(GEOquery)
  library(limma)
  library(ggplot2)
  library(data.table)
})

safe_mkdir <- function(d) if (!dir.exists(d)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

## ---- helpers ----
save_rds <- function(x, path) {
  safe_mkdir(dirname(path))
  saveRDS(x, path)
  msg('Saved: ', path)
}

write_csv <- function(df, path) {
  safe_mkdir(dirname(path))
  data.table::fwrite(df, path)
  msg('Saved: ', path)
}

## Extract expression matrix from ExpressionSet
get_expr_mat <- function(eset) {
  x <- Biobase::exprs(eset)
  as.matrix(x)
}

## If platform is microarray, log2 transform if needed
auto_log2 <- function(x) {
  qx <- quantile(x, probs = c(0, 0.25, 0.5, 0.75, 0.99, 1.0), na.rm = TRUE)
  # heuristic from GEO2R-like logic
  logc <- (qx[6] > 100) || (qx[5] > 50 && qx[2] > 0) || (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)
  if (logc) {
    msg('Applying log2(x + 1) transform (heuristic).')
    x <- log2(x + 1)
  } else {
    msg('No log2 transform applied (heuristic).')
  }
  x
}

## Map features to gene symbols if possible
## For many GEO series, feature data has 'Gene Symbol' or similar columns.
map_to_symbol <- function(expr, fdata) {
  col_candidates <- c('Gene Symbol','GENE_SYMBOL','gene_symbol','Symbol','SYMBOL','Gene.symbol','GENE','Gene Symbol(s)')
  hit <- intersect(col_candidates, colnames(fdata))
  if (length(hit) == 0) {
    # try fuzzy find
    idx <- grep('symbol|gene.*symbol', colnames(fdata), ignore.case = TRUE)
    if (length(idx) > 0) hit <- colnames(fdata)[idx[1]]
  }
  if (length(hit) == 0) {
    msg('No gene symbol column found in featureData; keeping probe-level matrix.')
    return(list(expr = expr, symbol_col = NA_character_, mapped = FALSE))
  }
  symcol <- hit[1]
  sym <- as.character(fdata[[symcol]])
  sym <- trimws(sym)
  sym[sym == '' | is.na(sym)] <- NA

  keep <- !is.na(sym)
  expr2 <- expr[keep, , drop = FALSE]
  sym2  <- sym[keep]

  # If multiple symbols separated by '///' or ';', keep first
  sym2 <- sub('\\s*///.*$', '', sym2)
  sym2 <- sub('\\s*;.*$', '', sym2)
  sym2 <- trimws(sym2)

  # Aggregate probes -> gene symbol by median
  dt <- data.table::as.data.table(expr2, keep.rownames = 'feature')
  dt[, symbol := sym2]
  dt <- dt[!is.na(symbol) & symbol != '']
  expr_cols <- setdiff(colnames(dt), c('feature','symbol'))
  dt_m <- dt[, lapply(.SD, median, na.rm = TRUE), by = symbol, .SDcols = expr_cols]
  mat <- as.matrix(dt_m[, ..expr_cols])
  rownames(mat) <- dt_m$symbol

  msg('Mapped to gene symbols using column: ', symcol, '; genes=', nrow(mat), '; samples=', ncol(mat))
  list(expr = mat, symbol_col = symcol, mapped = TRUE)
}

qc_plots <- function(expr, out_prefix) {
  safe_mkdir(dir_results)
  # PCA on samples
  x <- t(scale(t(expr), center = TRUE, scale = TRUE))
  x <- x[apply(x, 2, function(v) all(is.finite(v))), , drop = FALSE]
  pca <- prcomp(t(x), center = TRUE, scale. = FALSE)
  df <- data.frame(PC1 = pca$x[,1], PC2 = pca$x[,2], sample = rownames(pca$x))
  p <- ggplot(df, aes(PC1, PC2)) + geom_point() + theme_bw() + ggtitle('PCA (samples)')
  ggsave(filename = paste0(out_prefix, '_PCA.png'), plot = p, width = 6, height = 4, dpi = 200)
  msg('Saved QC PCA: ', paste0(out_prefix, '_PCA.png'))

  # Boxplot distribution (first 50 samples if huge)
  ns <- ncol(expr)
  pick <- if (ns > 50) 1:50 else 1:ns
  png(paste0(out_prefix, '_boxplot.png'), width = 1200, height = 700, res = 150)
  boxplot(expr[, pick, drop = FALSE], outline = FALSE, las = 2, cex.axis = 0.6, main = 'Expression distribution (subset samples)')
  dev.off()
  msg('Saved QC boxplot: ', paste0(out_prefix, '_boxplot.png'))
}

sepSMART_check <- function(expr, out_prefix) {
  have <- intersect(genes_sepSMART, rownames(expr))
  miss <- setdiff(genes_sepSMART, rownames(expr))
  df <- data.frame(gene = genes_sepSMART, present = genes_sepSMART %in% rownames(expr))
  write_csv(df, paste0(out_prefix, '_SepSMART_gene_presence.csv'))
  msg('SepSMART present: ', paste(have, collapse = ', '))
  if (length(miss) > 0) msg('SepSMART missing: ', paste(miss, collapse = ', '))

  if (length(have) >= 2) {
    sub <- expr[have, , drop = FALSE]
    # scale genes for heatmap-like visualization via base image
    z <- t(scale(t(sub)))
    z[!is.finite(z)] <- 0
    png(paste0(out_prefix, '_SepSMART_heatmap.png'), width = 900, height = 450, res = 150)
    op <- par(mar = c(5,6,3,1))
    image(t(z[nrow(z):1, , drop = FALSE]), axes = FALSE, main = 'SepSMART (z-score per gene)')
    axis(2, at = seq(0,1,length.out=nrow(z)), labels = rev(rownames(z)), las = 2, cex.axis = 0.8)
    par(op)
    dev.off()
    msg('Saved SepSMART heatmap: ', paste0(out_prefix, '_SepSMART_heatmap.png'))
  }
}

run_geo_pipeline <- function(gse_id, gpl_pick = NULL) {
  msg('=== GEO pipeline start: ', gse_id, ' ===')
  safe_mkdir(dir_raw); safe_mkdir(dir_clean); safe_mkdir(dir_results)

  # Download (getGEO may return a list if multiple GPLs)
  gse <- GEOquery::getGEO(gse_id, GSEMatrix = TRUE, getGPL = TRUE)

  if (is.list(gse) && length(gse) > 1) {
    msg('Multiple platforms detected: ', paste(names(gse), collapse = ', '))
    if (is.null(gpl_pick)) {
      msg('No gpl_pick provided; selecting the first platform: ', names(gse)[1])
      eset <- gse[[1]]
      picked <- names(gse)[1]
    } else {
      stopifnot(gpl_pick %in% names(gse))
      eset <- gse[[gpl_pick]]
      picked <- gpl_pick
    }
  } else {
    eset <- if (is.list(gse)) gse[[1]] else gse
    picked <- 'single_platform'
  }

  # Save raw ExpressionSet
  raw_path <- file.path(dir_raw, paste0(gse_id, '_', picked, '_ExpressionSet.rds'))
  save_rds(eset, raw_path)

  # Expression + pheno + feature data
  expr <- get_expr_mat(eset)
  pheno <- Biobase::pData(eset)
  fdata <- Biobase::fData(eset)

  # Basic cleaning
  expr <- auto_log2(expr)
  expr <- expr[apply(expr, 1, function(v) all(is.finite(v))), , drop = FALSE]

  # Map to gene symbols when possible
  mapped <- map_to_symbol(expr, fdata)
  expr2 <- mapped$expr

  # Save clean objects
  clean_expr_path <- file.path(dir_clean, paste0(gse_id, '_', picked, '_expr_clean.rds'))
  clean_pheno_path <- file.path(dir_clean, paste0(gse_id, '_', picked, '_pheno.rds'))
  save_rds(expr2, clean_expr_path)
  save_rds(pheno, clean_pheno_path)

  # Also save pheno as CSV for human inspection
  write_csv(as.data.frame(pheno), file.path(dir_clean, paste0(gse_id, '_', picked, '_pheno.csv')))

  # Output prefixes
  out_prefix <- file.path(dir_results, paste0('F0_', gse_id, '_', picked))

  # QC plots + SepSMART check
  qc_plots(expr2, out_prefix)
  sepSMART_check(expr2, out_prefix)

  msg('=== GEO pipeline done: ', gse_id, ' ===')
  invisible(list(expr = expr2, pheno = pheno, fdata = fdata, picked = picked))
}

msg('Loaded 02_geo_download_clean.R. Call run_geo_pipeline(GSE_ID).')
