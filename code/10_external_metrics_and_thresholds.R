## =========================================================
## 10_external_metrics_and_thresholds.R
## Purpose:
##   External metrics at an explicit threshold (Youden):
##   AUC + CI, threshold, Se/Sp/PPV/NPV/Acc
## Inputs:
##   02_data_clean/GSE26378_pheno_labeled.csv
##   02_data_clean/GSE28750_pheno_labeled.csv
##   02_data_clean/GSE26378_SepSMART_expr.csv (6x103, rows follow SEPSMART_GENES)
##   02_data_clean/GSE28750_SepSMART_expr.csv (6x41,  rows follow SEPSMART_GENES)
## Outputs:
##   04_results/T13_external_threshold_metrics.csv
## =========================================================

source('03_scripts/00_paths.R')
options(stringsAsFactors = FALSE)

SEPSMART_GENES <- c('RETN','MCEMP1','CYP1B1','S100A12','S100A8','HK3')

pkgs <- c('readr','dplyr','stringr','tibble','pROC')
missing <- pkgs[!vapply(pkgs, requireNamespace, quietly=TRUE, FUN.VALUE=logical(1))]
if (length(missing)>0) {
  if (requireNamespace('renv', quietly=TRUE)) try(renv::install(missing), silent=TRUE)
  missing2 <- missing[!vapply(missing, requireNamespace, quietly=TRUE, FUN.VALUE=logical(1))]
  if (length(missing2)>0) install.packages(missing2)
}
suppressPackageStartupMessages({
  library(readr); library(dplyr); library(stringr); library(tibble); library(pROC)
})

# ---- robust reader: your expr csv has NO gene column, rows follow SEPSMART_GENES order ----
read_expr_with_genes <- function(path, genes=SEPSMART_GENES) {
  stopifnot(file.exists(path))
  x <- readr::read_csv(path, show_col_types=FALSE)
  # if gene column exists, use it; otherwise assume row order = genes
  if (ncol(x) >= 2 && any(tolower(names(x)) %in% c('gene','symbol'))) {
    gene_col <- names(x)[tolower(names(x)) %in% c('gene','symbol')][1]
    mat <- as.matrix(x[, setdiff(names(x), gene_col), drop=FALSE])
    rownames(mat) <- toupper(as.character(x[[gene_col]]))
  } else {
    mat <- as.matrix(x)
    if (nrow(mat) != length(genes)) {
      stop('Expr file rows != length(SEPSMART_GENES). Please re-export with gene names. File: ', path)
    }
    rownames(mat) <- genes
  }
  storage.mode(mat) <- 'numeric'
  mat
}

# ---- score: gene-wise z within cohort then mean ----
score_sepsmart <- function(gene_mat, genes=SEPSMART_GENES) {
  g <- intersect(genes, rownames(gene_mat))
  if (length(g) < 5) stop('Too few genes present (<5).')
  m <- gene_mat[g, , drop=FALSE]
  mz <- t(scale(t(m)))
  mz[is.na(mz)] <- 0
  as.numeric(colMeans(mz))
}

confusion_metrics <- function(y, score, thr) {
  pred <- ifelse(score >= thr, 1L, 0L)
  tp <- sum(pred==1 & y==1); tn <- sum(pred==0 & y==0)
  fp <- sum(pred==1 & y==0); fn <- sum(pred==0 & y==1)
  se <- ifelse((tp+fn)>0, tp/(tp+fn), NA_real_)
  sp <- ifelse((tn+fp)>0, tn/(tn+fp), NA_real_)
  ppv <- ifelse((tp+fp)>0, tp/(tp+fp), NA_real_)
  npv <- ifelse((tn+fn)>0, tn/(tn+fn), NA_real_)
  acc <- (tp+tn)/length(y)
  tibble(tp=tp, tn=tn, fp=fp, fn=fn, sensitivity=se, specificity=sp, ppv=ppv, npv=npv, accuracy=acc)
}

eval_endpoint <- function(endpoint_name, y, score) {
  roc <- pROC::roc(y, score, quiet=TRUE)
  auc <- as.numeric(pROC::auc(roc))
  ci <- as.numeric(pROC::ci.auc(roc))
  # Youden threshold
  thr <- as.numeric(pROC::coords(roc, x='best', best.method='youden', ret='threshold', transpose=FALSE))
  met <- confusion_metrics(y, score, thr)
  tibble(endpoint=endpoint_name, n=length(y), pos=sum(y==1), neg=sum(y==0),
         auc=auc, auc_ci_low=ci[1], auc_ci_high=ci[3], threshold_youden=thr) %>%
    bind_cols(met)
}

# ---- Load labeled pheno ----
p26378 <- file.path(dir_clean, 'GSE26378_pheno_labeled.csv')
p28750 <- file.path(dir_clean, 'GSE28750_pheno_labeled.csv')
stopifnot(file.exists(p26378), file.exists(p28750))
ph26378 <- readr::read_csv(p26378, show_col_types=FALSE) %>% mutate(.sid=as.character(geo_accession))
ph28750 <- readr::read_csv(p28750, show_col_types=FALSE) %>% mutate(.sid=as.character(geo_accession))

# ---- Load expr ----
e26378 <- read_expr_with_genes(file.path(dir_clean,'GSE26378_SepSMART_expr.csv'))
e28750 <- read_expr_with_genes(file.path(dir_clean,'GSE28750_SepSMART_expr.csv'))

# ---- Align + score ----
common26378 <- intersect(colnames(e26378), ph26378$.sid)
common28750 <- intersect(colnames(e28750), ph28750$.sid)
stopifnot(length(common26378)>0, length(common28750)>0)

ph26378 <- ph26378 %>% filter(.sid %in% common26378) %>% arrange(match(.sid, common26378))
ph28750 <- ph28750 %>% filter(.sid %in% common28750) %>% arrange(match(.sid, common28750))

s26378 <- score_sepsmart(e26378[, common26378, drop=FALSE])
s28750 <- score_sepsmart(e28750[, common28750, drop=FALSE])

# ---- Endpoints ----
y26378_dx <- ifelse(ph26378$dx_group=='Septic_shock', 1L, 0L)
out <- list()
out[[1]] <- eval_endpoint('GSE26378 Diagnosis (Septic_shock vs Control)', y26378_dx, s26378)

if ('mort' %in% names(ph26378) && any(!is.na(ph26378$mort))) {
  keep <- !is.na(ph26378$mort)
  out[[2]] <- eval_endpoint('GSE26378 Mortality (Nonsurvivor vs Survivor)', as.integer(ph26378$mort[keep]), s26378[keep])
}

keepA <- !is.na(ph28750$dx_A)
yA <- ifelse(ph28750$dx_A[keepA]=='Sepsis', 1L, 0L)
out[[3]] <- eval_endpoint('GSE28750(A) Diagnosis (Sepsis vs Post_surgical)', yA, s28750[keepA])

keepB <- !is.na(ph28750$dx_B)
yB <- ifelse(ph28750$dx_B[keepB]=='Sepsis', 1L, 0L)
out[[4]] <- eval_endpoint('GSE28750(B) Diagnosis (Sepsis vs All controls)', yB, s28750[keepB])

T13 <- bind_rows(out)
out_path <- file.path(dir_results, 'T13_external_threshold_metrics.csv')
readr::write_csv(T13, out_path)
cat('✅ Saved: ', out_path, '\n', sep='')
print(T13)
