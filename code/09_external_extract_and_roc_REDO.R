## =========================================================
## 09_external_extract_and_roc_REDO.R
## Rebuild external 6-gene expression exports and ROC summary
## Key fixes:
##   1) auto_log2() before collapse/export
##   2) harmonize alias C19orf59 -> MCEMP1
##   3) enforce locked gene order in exported expression CSV
##   4) sample-level score_sepsmart()
## =========================================================

source(file.path("/Users/felix/Documents/SepSMART_restart/03_scripts", "00_paths.R"))

suppressPackageStartupMessages({
  library(Biobase)
  library(readr)
  library(dplyr)
  library(ggplot2)
  library(pROC)
  library(data.table)
})

options(stringsAsFactors = FALSE, scipen = 999)

if (!exists("project_root")) {
  project_root <- "/Users/felix/Documents/SepSMART_restart"
}
if (!exists("dir_clean")) {
  dir_clean <- file.path(project_root, "02_data_clean")
}
if (!exists("dir_results")) {
  dir_results <- file.path(project_root, "04_results")
}
if (!exists("dir_raw")) {
  dir_raw <- file.path(project_root, "01_data_raw")
}

dir.create(dir_clean, showWarnings = FALSE, recursive = TRUE)
dir.create(dir_results, showWarnings = FALSE, recursive = TRUE)

msg <- function(...) cat(format(Sys.time(), "[%Y-%m-%d %H:%M:%S] "), ..., "\n")

SEPSMART_GENES <- c("RETN", "MCEMP1", "CYP1B1", "S100A12", "S100A8", "HK3")

## -----------------------------
## helpers
## -----------------------------
auto_log2 <- function(x) {
  qx <- quantile(x, probs = c(0, 0.25, 0.5, 0.75, 0.99, 1.0), na.rm = TRUE)
  logc <- (qx[6] > 100) || (qx[5] > 50 && qx[2] > 0) || (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)
  if (logc) {
    msg("Applying log2(x + 1) transform (heuristic).")
    x <- log2(x + 1)
  } else {
    msg("No log2 transform applied (heuristic).")
  }
  x
}

probe_symbol_map <- function(eset) {
  fdata <- Biobase::fData(eset)
  
  col_candidates <- c(
    "Gene Symbol", "GENE_SYMBOL", "gene_symbol", "Symbol", "SYMBOL",
    "Gene.symbol", "GENE", "Gene Symbol(s)"
  )
  hit <- intersect(col_candidates, colnames(fdata))
  if (length(hit) == 0) {
    idx <- grep("symbol|gene.*symbol", colnames(fdata), ignore.case = TRUE)
    if (length(idx) > 0) hit <- colnames(fdata)[idx[1]]
  }
  
  if (length(hit) == 0) {
    msg("No gene symbol column found in featureData; keeping feature names as-is.")
    sym <- rownames(fdata)
    note <- "No symbol column; used feature names"
  } else {
    symcol <- hit[1]
    sym <- as.character(fdata[[symcol]])
    note <- paste0("fData:", symcol)
  }
  
  sym <- trimws(sym)
  sym[sym == "" | is.na(sym)] <- NA_character_
  
  map <- data.frame(
    probe = rownames(fdata),
    symbol = sym,
    stringsAsFactors = FALSE
  )
  
  list(map = map, note = note)
}

collapse_to_gene <- function(expr, map_df) {
  stopifnot(is.matrix(expr))
  stopifnot(all(rownames(expr) %in% map_df$probe))
  
  ord <- match(rownames(expr), map_df$probe)
  map_df2 <- map_df[ord, , drop = FALSE]
  
  keep <- !is.na(map_df2$symbol) & nzchar(map_df2$symbol)
  expr2 <- expr[keep, , drop = FALSE]
  sym2 <- toupper(map_df2$symbol[keep])
  
  dt <- data.table::as.data.table(expr2, keep.rownames = "probe")
  dt[, symbol := sym2]
  expr_cols <- setdiff(colnames(dt), c("probe", "symbol"))
  
  dt_m <- dt[, lapply(.SD, median, na.rm = TRUE), by = symbol, .SDcols = expr_cols]
  mat <- as.matrix(dt_m[, ..expr_cols])
  rownames(mat) <- dt_m$symbol
  mode(mat) <- "numeric"
  mat
}

harmonize_aliases <- function(gene_mat) {
  rn <- rownames(gene_mat)
  rn[rn == "C19ORF59"] <- "MCEMP1"
  rownames(gene_mat) <- rn
  gene_mat
}

prepare_out_expr <- function(gene_mat, genes = SEPSMART_GENES) {
  gene_mat <- harmonize_aliases(gene_mat)
  miss <- setdiff(genes, rownames(gene_mat))
  if (length(miss) > 0) {
    stop("Missing required genes after collapse/harmonization: ", paste(miss, collapse = ", "))
  }
  gene_mat[genes, , drop = FALSE]
}

## sample-level exploratory SepSMART score
## Note: this is not the locked GLM deployment score.
score_sepsmart <- function(gene_mat) {
  gene_mat <- harmonize_aliases(gene_mat)
  miss <- setdiff(SEPSMART_GENES, rownames(gene_mat))
  if (length(miss) > 0) {
    stop("Missing genes for score_sepsmart: ", paste(miss, collapse = ", "))
  }
  
  X <- gene_mat[SEPSMART_GENES, , drop = FALSE]   # genes x samples
  Z <- t(scale(t(X)))
  Z[is.na(Z)] <- 0
  
  w <- c(
    RETN    =  1,
    MCEMP1  = -1,
    CYP1B1  = -1,
    S100A12 =  1,
    S100A8  = -1,
    HK3     =  1
  )
  
  ## t(Z): samples x genes
  score <- as.numeric(rowMeans(t(Z) * w[colnames(t(Z))]))
  score
}

eval_roc <- function(y, score) {
  roc_obj <- pROC::roc(y, score, quiet = TRUE, direction = "<")
  ci <- as.numeric(pROC::ci.auc(roc_obj))
  list(
    roc = roc_obj,
    auc = as.numeric(pROC::auc(roc_obj)),
    ci_low = ci[1],
    ci_high = ci[3]
  )
}

plot_roc <- function(roc_obj, title, out_png) {
  png(out_png, width = 900, height = 700, res = 150)
  plot(roc_obj, main = title, legacy.axes = TRUE)
  abline(a = 0, b = 1, lty = 2)
  dev.off()
  message("Saved: ", out_png)
}

write_external_expr_csv <- function(out_expr, out_path, genes = SEPSMART_GENES) {
  tmp <- as.data.frame(out_expr)
  if (nrow(tmp) != length(genes)) {
    stop("out_expr rows != length(SEPSMART_GENES)")
  }
  tmp2 <- cbind(gene = genes, tmp)
  readr::write_csv(tmp2, out_path)
  msg("Saved: ", out_path)
}

## -----------------------------
## cohort definitions
## -----------------------------
COHORTS <- list(
  GSE26378 = list(
    pheno_labeled = file.path(dir_clean, "GSE26378_pheno_labeled.csv"),
    raw_eset = file.path(dir_raw, "GSE26378", "GSE26378_eset.rds"),
    dx_field = "dx_group",
    dx_pos = "Septic_shock",
    mort_field = "mort"
  ),
  GSE28750 = list(
    pheno_labeled = file.path(dir_clean, "GSE28750_pheno_labeled.csv"),
    raw_eset = file.path(dir_raw, "GSE28750", "GSE28750_eset.rds")
  )
)

summary_rows <- list()

## -----------------------------
## GSE26378
## -----------------------------
{
  info <- COHORTS$GSE26378
  ph <- readr::read_csv(info$pheno_labeled, show_col_types = FALSE)
  eset <- readRDS(info$raw_eset)
  
  expr <- Biobase::exprs(eset)
  expr <- auto_log2(expr)
  
  mp <- probe_symbol_map(eset)
  gene_mat <- collapse_to_gene(expr, mp$map)
  gene_mat <- harmonize_aliases(gene_mat)
  
  genes_have <- intersect(SEPSMART_GENES, rownames(gene_mat))
  out_expr <- prepare_out_expr(gene_mat, SEPSMART_GENES)
  
  out_path <- file.path(dir_clean, "GSE26378_SepSMART_expr.csv")
  write_external_expr_csv(out_expr, out_path, SEPSMART_GENES)
  
  samp <- colnames(gene_mat)
  ph <- ph %>% mutate(.sample_id = as.character(geo_accession))
  common <- intersect(samp, ph$.sample_id)
  ph2 <- ph %>% filter(.sample_id %in% common) %>% arrange(match(.sample_id, common))
  
  score <- score_sepsmart(gene_mat[, common, drop = FALSE])
  y_dx <- ifelse(ph2[[info$dx_field]] == info$dx_pos, 1L, 0L)
  
  cat("GSE26378: length(common) =", length(common), "\n")
  cat("GSE26378: nrow(ph2) =", nrow(ph2), "\n")
  cat("GSE26378: length(score) =", length(score), "\n")
  cat("GSE26378: length(y_dx) =", length(y_dx), "\n")
  
  rdx <- eval_roc(y_dx, score)
  plot_roc(
    rdx$roc,
    sprintf("GSE26378 Diagnosis ROC (AUC=%.3f)", rdx$auc),
    file.path(dir_results, "F8_GSE26378_dx_ROC.png")
  )
  
  summary_rows[[length(summary_rows) + 1]] <- tibble(
    cohort = "GSE26378",
    endpoint = "Diagnosis (Septic_shock vs Control)",
    n = length(y_dx),
    pos = sum(y_dx == 1),
    neg = sum(y_dx == 0),
    genes_present = length(genes_have),
    genes_missing = paste(setdiff(SEPSMART_GENES, genes_have), collapse = ";"),
    auc = rdx$auc,
    auc_ci_low = rdx$ci_low,
    auc_ci_high = rdx$ci_high,
    mapping_note = mp$note
  )
  
  if ("mort" %in% colnames(ph2) && any(!is.na(ph2$mort))) {
    keep <- !is.na(ph2$mort)
    y_m <- as.integer(ph2$mort[keep])
    rmo <- eval_roc(y_m, score[keep])
    plot_roc(
      rmo$roc,
      sprintf("GSE26378 Mortality ROC (AUC=%.3f)", rmo$auc),
      file.path(dir_results, "F9_GSE26378_mortality_ROC.png")
    )
    
    summary_rows[[length(summary_rows) + 1]] <- tibble(
      cohort = "GSE26378",
      endpoint = "Mortality (Nonsurvivor vs Survivor)",
      n = length(y_m),
      pos = sum(y_m == 1),
      neg = sum(y_m == 0),
      genes_present = length(genes_have),
      genes_missing = paste(setdiff(SEPSMART_GENES, genes_have), collapse = ";"),
      auc = rmo$auc,
      auc_ci_low = rmo$ci_low,
      auc_ci_high = rmo$ci_high,
      mapping_note = mp$note
    )
  }
}

## -----------------------------
## GSE28750
## -----------------------------
{
  info <- COHORTS$GSE28750
  ph <- readr::read_csv(info$pheno_labeled, show_col_types = FALSE)
  eset <- readRDS(info$raw_eset)
  
  expr <- Biobase::exprs(eset)
  expr <- auto_log2(expr)
  
  mp <- probe_symbol_map(eset)
  gene_mat <- collapse_to_gene(expr, mp$map)
  gene_mat <- harmonize_aliases(gene_mat)
  
  genes_have <- intersect(SEPSMART_GENES, rownames(gene_mat))
  out_expr <- prepare_out_expr(gene_mat, SEPSMART_GENES)
  
  out_path <- file.path(dir_clean, "GSE28750_SepSMART_expr.csv")
  write_external_expr_csv(out_expr, out_path, SEPSMART_GENES)
  
  ph <- ph %>% mutate(.sample_id = as.character(geo_accession))
  samp <- colnames(gene_mat)
  common <- intersect(samp, ph$.sample_id)
  ph2 <- ph %>% filter(.sample_id %in% common) %>% arrange(match(.sample_id, common))
  
  score <- score_sepsmart(gene_mat[, common, drop = FALSE])
  
  cat("GSE28750: length(common) =", length(common), "\n")
  cat("GSE28750: nrow(ph2) =", nrow(ph2), "\n")
  cat("GSE28750: length(score) =", length(score), "\n")
  
  if ("dx_A" %in% colnames(ph2)) {
    keepA <- !is.na(ph2$dx_A)
    yA <- ifelse(ph2$dx_A[keepA] == "Sepsis", 1L, 0L)
    rA <- eval_roc(yA, score[keepA])
    plot_roc(
      rA$roc,
      sprintf("GSE28750A Diagnosis ROC (AUC=%.3f)", rA$auc),
      file.path(dir_results, "F8_GSE28750A_dx_ROC.png")
    )
    
    summary_rows[[length(summary_rows) + 1]] <- tibble(
      cohort = "GSE28750",
      endpoint = "Diagnosis A (Sepsis vs Post_surgical)",
      n = length(yA),
      pos = sum(yA == 1),
      neg = sum(yA == 0),
      genes_present = length(genes_have),
      genes_missing = paste(setdiff(SEPSMART_GENES, genes_have), collapse = ";"),
      auc = rA$auc,
      auc_ci_low = rA$ci_low,
      auc_ci_high = rA$ci_high,
      mapping_note = mp$note
    )
  }
  
  if ("dx_B" %in% colnames(ph2)) {
    keepB <- !is.na(ph2$dx_B)
    yB <- ifelse(ph2$dx_B[keepB] == "Sepsis", 1L, 0L)
    rB <- eval_roc(yB, score[keepB])
    plot_roc(
      rB$roc,
      sprintf("GSE28750B Diagnosis ROC (AUC=%.3f)", rB$auc),
      file.path(dir_results, "F8_GSE28750B_dx_ROC.png")
    )
    
    summary_rows[[length(summary_rows) + 1]] <- tibble(
      cohort = "GSE28750",
      endpoint = "Diagnosis B (Sepsis vs All controls)",
      n = length(yB),
      pos = sum(yB == 1),
      neg = sum(yB == 0),
      genes_present = length(genes_have),
      genes_missing = paste(setdiff(SEPSMART_GENES, genes_have), collapse = ";"),
      auc = rB$auc,
      auc_ci_low = rB$ci_low,
      auc_ci_high = rB$ci_high,
      mapping_note = mp$note
    )
  }
}

## -----------------------------
## save summary
## -----------------------------
T12 <- bind_rows(summary_rows)
readr::write_csv(T12, file.path(dir_results, "T12_external_ROC_summary.csv"))
cat("✅ Done. Summary: ", file.path(dir_results, "T12_external_ROC_summary.csv"), "\n", sep = "")
print(T12)