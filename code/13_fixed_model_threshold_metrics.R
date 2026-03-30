## 13_fixed_model_threshold_metrics.R
## Fixed-threshold external metrics for locked GLM
## Revised:
##   1) exact binomial 95% CI for sensitivity/specificity/PPV/NPV/accuracy
##   2) exports original T17 and new T17b with CIs

source(file.path("/Users/felix/Documents/SepSMART_restart/03_scripts", "00_paths.R"))

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tibble)
  library(pROC)
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

SEPSMART_GENES <- c("RETN","MCEMP1","CYP1B1","S100A12","S100A8","HK3")

read_expr_with_genes <- function(path, genes = SEPSMART_GENES) {
  stopifnot(file.exists(path))
  x <- readr::read_csv(path, show_col_types = FALSE)
  if (ncol(x) >= 2 && any(tolower(names(x)) %in% c("gene","symbol"))) {
    gene_col <- names(x)[tolower(names(x)) %in% c("gene","symbol")][1]
    mat <- as.matrix(x[, setdiff(names(x), gene_col), drop = FALSE])
    rownames(mat) <- toupper(as.character(x[[gene_col]]))
  } else {
    mat <- as.matrix(x)
    if (nrow(mat) != length(genes)) stop("Expr rows != 6 and no Gene column. File: ", path)
    rownames(mat) <- genes
  }
  storage.mode(mat) <- "numeric"
  mat
}

apply_locked_model <- function(expr_mat, scaling_df, coef_tbl) {
  scaling_df <- scaling_df %>% mutate(gene = toupper(gene))
  genes_use <- scaling_df$gene
  
  miss <- setdiff(genes_use, rownames(expr_mat))
  if (length(miss) > 0) stop("External expr missing gene(s) required by scaling: ", paste(miss, collapse = ";"))
  
  X <- t(expr_mat[genes_use, , drop = FALSE])  # samples x genes
  
  for (j in seq_along(genes_use)) {
    g <- genes_use[j]
    X[, g] <- (X[, g] - scaling_df$mean[j]) / scaling_df$sd[j]
  }
  
  b0 <- coef_tbl$coefficient[coef_tbl$term == "(Intercept)"][1]
  b  <- coef_tbl$coefficient[match(genes_use, coef_tbl$term)]
  if (any(is.na(b))) stop("Coefficient table missing gene terms: ", paste(genes_use[is.na(b)], collapse = ";"))
  
  lp <- as.numeric(b0 + X %*% b)
  p  <- 1 / (1 + exp(-lp))
  list(lp = lp, prob = p)
}

binom_ci <- function(x, n, conf.level = 0.95) {
  if (is.na(x) || is.na(n) || n <= 0) {
    return(c(estimate = NA_real_, low = NA_real_, high = NA_real_))
  }
  bt <- binom.test(x, n, conf.level = conf.level)
  c(
    estimate = unname(bt$estimate),
    low = bt$conf.int[1],
    high = bt$conf.int[2]
  )
}

safe_div <- function(num, den) {
  ifelse(den > 0, num / den, NA_real_)
}

conf_mat_metrics <- function(y, p, thr) {
  pred <- ifelse(p >= thr, 1L, 0L)
  
  tp <- sum(pred == 1 & y == 1)
  tn <- sum(pred == 0 & y == 0)
  fp <- sum(pred == 1 & y == 0)
  fn <- sum(pred == 0 & y == 1)
  
  sens <- safe_div(tp, tp + fn)
  spec <- safe_div(tn, tn + fp)
  ppv  <- safe_div(tp, tp + fp)
  npv  <- safe_div(tn, tn + fn)
  acc  <- safe_div(tp + tn, length(y))
  
  sens_ci <- binom_ci(tp, tp + fn)
  spec_ci <- binom_ci(tn, tn + fp)
  ppv_ci  <- binom_ci(tp, tp + fp)
  npv_ci  <- binom_ci(tn, tn + fn)
  acc_ci  <- binom_ci(tp + tn, length(y))
  
  tibble(
    tp = tp, tn = tn, fp = fp, fn = fn,
    sensitivity = sens,
    sensitivity_ci_low = sens_ci["low"],
    sensitivity_ci_high = sens_ci["high"],
    specificity = spec,
    specificity_ci_low = spec_ci["low"],
    specificity_ci_high = spec_ci["high"],
    ppv = ppv,
    ppv_ci_low = ppv_ci["low"],
    ppv_ci_high = ppv_ci["high"],
    npv = npv,
    npv_ci_low = npv_ci["low"],
    npv_ci_high = npv_ci["high"],
    accuracy = acc,
    accuracy_ci_low = acc_ci["low"],
    accuracy_ci_high = acc_ci["high"]
  )
}

## ---- Locked artifacts ----
coef_path <- file.path(dir_results, "T14_GSE65682_dx_model_coefficients.csv")
sc_path   <- file.path(dir_results, "T15_GSE65682_dx_scaling_params.csv")
stopifnot(file.exists(coef_path), file.exists(sc_path))

coef_all <- readr::read_csv(coef_path, show_col_types = FALSE)
scaling  <- readr::read_csv(sc_path, show_col_types = FALSE)
coef_glm <- coef_all %>% filter(model == "GLM") %>% select(term, coefficient)

## ---- Training CV Youden threshold ----
tr <- readr::read_csv(file.path(dir_results, "T14b_GSE65682_dx_internal_cv_predictions.csv"), show_col_types = FALSE)
roc_tr <- pROC::roc(tr$y, tr$pred_glm_cv, quiet = TRUE)
thr_train_youden <- as.numeric(
  pROC::coords(
    roc_tr,
    x = "best",
    best.method = "youden",
    ret = "threshold",
    transpose = FALSE
  )
)
cat("[INFO] Training CV Youden threshold (GLM)=", thr_train_youden, "\n", sep = "")

## ---- External inputs ----
ph26378 <- readr::read_csv(file.path(dir_clean, "GSE26378_pheno_labeled.csv"), show_col_types = FALSE) %>%
  mutate(.sid = as.character(geo_accession))
ph28750 <- readr::read_csv(file.path(dir_clean, "GSE28750_pheno_labeled.csv"), show_col_types = FALSE) %>%
  mutate(.sid = as.character(geo_accession))

e26378 <- read_expr_with_genes(file.path(dir_clean, "GSE26378_SepSMART_expr.csv"))
e28750 <- read_expr_with_genes(file.path(dir_clean, "GSE28750_SepSMART_expr.csv"))

c1 <- intersect(colnames(e26378), ph26378$.sid)
c2 <- intersect(colnames(e28750), ph28750$.sid)

ph26378 <- ph26378 %>% filter(.sid %in% c1) %>% arrange(match(.sid, c1))
ph28750 <- ph28750 %>% filter(.sid %in% c2) %>% arrange(match(.sid, c2))

e26378 <- e26378[, c1, drop = FALSE]
e28750 <- e28750[, c2, drop = FALSE]

## ---- Endpoints ----
y26378 <- ifelse(ph26378$dx_group == "Septic_shock", 1L, 0L)

keepA <- !is.na(ph28750$dx_A)
yA <- ifelse(ph28750$dx_A[keepA] == "Sepsis", 1L, 0L)

keepB <- !is.na(ph28750$dx_B)
yB <- ifelse(ph28750$dx_B[keepB] == "Sepsis", 1L, 0L)

## ---- Locked GLM predictions ----
pred26378 <- apply_locked_model(e26378, scaling, coef_glm)$prob
pred28750 <- apply_locked_model(e28750, scaling, coef_glm)$prob

## ---- Build T17 / T17b ----
rows <- list()

add_thr <- function(model, cohort, endpoint, y, p, thr, thr_type) {
  met <- conf_mat_metrics(y, p, thr)
  
  tibble(
    model = model,
    cohort = cohort,
    endpoint = endpoint,
    threshold_type = thr_type,
    threshold = thr,
    n = length(y),
    pos = sum(y == 1),
    neg = sum(y == 0)
  ) %>%
    bind_cols(met)
}

rows[[1]] <- add_thr("GLM", "GSE26378", "Diagnosis (Septic_shock vs Control)", y26378, pred26378, 0.5, "Fixed_0.5")
rows[[2]] <- add_thr("GLM", "GSE26378", "Diagnosis (Septic_shock vs Control)", y26378, pred26378, thr_train_youden, "TrainCV_Youden")

rows[[3]] <- add_thr("GLM", "GSE28750", "Diagnosis A (Sepsis vs Post_surgical)", yA, pred28750[keepA], 0.5, "Fixed_0.5")
rows[[4]] <- add_thr("GLM", "GSE28750", "Diagnosis A (Sepsis vs Post_surgical)", yA, pred28750[keepA], thr_train_youden, "TrainCV_Youden")

rows[[5]] <- add_thr("GLM", "GSE28750", "Diagnosis B (Sepsis vs All controls)", yB, pred28750[keepB], 0.5, "Fixed_0.5")
rows[[6]] <- add_thr("GLM", "GSE28750", "Diagnosis B (Sepsis vs All controls)", yB, pred28750[keepB], thr_train_youden, "TrainCV_Youden")

T17b <- bind_rows(rows)

## backward-compatible T17 (point estimates only)
T17 <- T17b %>%
  select(
    model, cohort, endpoint, threshold_type, threshold,
    n, pos, neg, tp, tn, fp, fn,
    sensitivity, specificity, ppv, npv, accuracy
  )

out_t17 <- file.path(dir_results, "T17_fixed_model_external_threshold_metrics.csv")
out_t17b <- file.path(dir_results, "T17b_fixed_model_external_threshold_metrics_with_CI.csv")

readr::write_csv(T17, out_t17)
readr::write_csv(T17b, out_t17b)

cat("✅ Saved: ", out_t17, "\n", sep = "")
cat("✅ Saved: ", out_t17b, "\n", sep = "")

print(T17b)