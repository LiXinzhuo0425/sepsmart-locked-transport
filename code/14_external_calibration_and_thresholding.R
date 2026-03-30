## 14_external_calibration_and_thresholding.R
## External calibration / thresholding for locked SepSMART GLM
## Aligned to 12_validate_fixed_model_external.R scoring logic
## Revised:
##   1) identical locked-model scoring logic as T16 script
##   2) standard calibration intercept/slope definition
##   3) epsilon documented as 1e-6

source(file.path("/Users/felix/Documents/SepSMART_restart/03_scripts", "00_paths.R"))

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tibble)
  library(pROC)
  library(ggplot2)
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
  if (length(miss) > 0) {
    stop("External expr missing gene(s) required by scaling: ", paste(miss, collapse = ";"))
  }
  
  X <- t(expr_mat[genes_use, , drop = FALSE])  # samples x genes
  
  for (j in seq_along(genes_use)) {
    g <- genes_use[j]
    X[, g] <- (X[, g] - scaling_df$mean[j]) / scaling_df$sd[j]
  }
  
  b0 <- coef_tbl$coefficient[coef_tbl$term == "(Intercept)"][1]
  b  <- coef_tbl$coefficient[match(genes_use, coef_tbl$term)]
  if (any(is.na(b))) {
    stop("Coefficient table missing gene terms: ", paste(genes_use[is.na(b)], collapse = ";"))
  }
  
  lp <- as.numeric(b0 + X %*% b)
  p  <- 1 / (1 + exp(-lp))
  list(lp = lp, prob = p)
}

brier <- function(y, p) {
  mean((as.numeric(p) - as.numeric(y))^2, na.rm = TRUE)
}

calibration_fit <- function(y, p) {
  eps <- 1e-6
  p2 <- pmin(pmax(as.numeric(p), eps), 1 - eps)
  logit <- log(p2 / (1 - p2))
  
  dat <- data.frame(y = as.numeric(y), logit = logit)
  dat <- dat[complete.cases(dat), , drop = FALSE]
  
  fit_int <- glm(y ~ 1 + offset(logit), family = binomial(), data = dat)
  co_int <- coef(summary(fit_int))
  
  fit_slo <- glm(y ~ logit, family = binomial(), data = dat)
  co_slo <- coef(summary(fit_slo))
  
  tibble(
    cal_intercept = unname(co_int[1, "Estimate"]),
    cal_slope     = unname(co_slo[2, "Estimate"]),
    intercept_se  = unname(co_int[1, "Std. Error"]),
    slope_se      = unname(co_slo[2, "Std. Error"])
  )
}

calibration_curve <- function(y, p, bins = 10) {
  tibble(y = y, p = p) %>%
    mutate(bin = ntile(p, bins)) %>%
    group_by(bin) %>%
    summarise(
      p_mean = mean(p, na.rm = TRUE),
      y_mean = mean(y, na.rm = TRUE),
      n = n(),
      .groups = "drop"
    )
}

plot_score_dist <- function(df, title, out_png) {
  p <- ggplot(df, aes(x = factor(y), y = prob)) +
    geom_boxplot(outlier.shape = NA, linewidth = 0.7) +
    geom_jitter(width = 0.15, height = 0, alpha = 0.5, size = 1) +
    theme_bw() +
    labs(
      title = title,
      x = "Class (0=Control, 1=Case)",
      y = "Predicted risk (locked GLM)"
    )
  ggsave(out_png, p, width = 5.8, height = 4.6, dpi = 220)
}

plot_calibration <- function(curve_df, title, out_png) {
  p <- ggplot(curve_df, aes(x = p_mean, y = y_mean, size = n)) +
    geom_point(alpha = 0.85) +
    geom_line(linewidth = 0.8) +
    geom_abline(slope = 1, intercept = 0, linetype = 2) +
    theme_bw() +
    labs(
      title = title,
      x = "Mean predicted risk",
      y = "Observed event rate",
      size = "Bin n"
    ) +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1))
  ggsave(out_png, p, width = 5.8, height = 4.8, dpi = 220)
}

find_thr_youden <- function(y, p) {
  roc <- pROC::roc(y, p, quiet = TRUE)
  as.numeric(
    pROC::coords(
      roc,
      x = "best",
      best.method = "youden",
      ret = "threshold",
      transpose = FALSE
    )
  )
}

find_thr_spec <- function(y, p, target_spec = 0.90) {
  roc <- pROC::roc(y, p, quiet = TRUE)
  df <- pROC::coords(
    roc,
    x = "all",
    ret = c("threshold", "specificity", "sensitivity"),
    transpose = FALSE
  ) %>%
    as.data.frame() %>%
    rename(
      thr = threshold,
      sp = specificity,
      se = sensitivity
    ) %>%
    filter(!is.na(thr), !is.na(sp), !is.na(se)) %>%
    filter(sp >= target_spec) %>%
    arrange(desc(se), desc(sp))
  if (nrow(df) == 0) return(NA_real_)
  df$thr[1]
}

## ---- Load locked artifacts ----
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
cat("[INFO] Train CV Youden threshold=", thr_train_youden, "\n", sep = "")

## ---- External files ----
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

## ---- Locked GLM predictions: EXACTLY aligned with T16 logic ----
pred26378 <- apply_locked_model(e26378, scaling, coef_glm)
pred28750 <- apply_locked_model(e28750, scaling, coef_glm)

p26378 <- pred26378$prob
p28750 <- pred28750$prob

cat("length(y26378) =", length(y26378), "\n")
cat("length(p26378) =", length(p26378), "\n")
cat("length(yA)     =", length(yA), "\n")
cat("length(p28750[keepA]) =", length(p28750[keepA]), "\n")
cat("length(yB)     =", length(yB), "\n")
cat("length(p28750[keepB]) =", length(p28750[keepB]), "\n")

## ---- Calibration summary ----
add_cal <- function(cohort, endpoint, y, p) {
  roc <- pROC::roc(y, p, quiet = TRUE)
  auc <- as.numeric(pROC::auc(roc))
  ci <- as.numeric(pROC::ci.auc(roc))
  cf <- calibration_fit(y, p)
  
  tibble(
    cohort = cohort,
    endpoint = endpoint,
    n = length(y),
    pos = sum(y == 1),
    neg = sum(y == 0),
    auc = auc,
    auc_ci_low = ci[1],
    auc_ci_high = ci[3],
    brier = brier(y, p)
  ) %>%
    bind_cols(cf)
}

T18 <- bind_rows(
  add_cal("GSE26378", "Diagnosis (Septic_shock vs Control)", y26378, p26378),
  add_cal("GSE28750", "Diagnosis A (Sepsis vs Post_surgical)", yA, p28750[keepA]),
  add_cal("GSE28750", "Diagnosis B (Sepsis vs All controls)", yB, p28750[keepB])
)

readr::write_csv(T18, file.path(dir_results, "T18_external_calibration_summary.csv"))
cat("✅ Saved: 04_results/T18_external_calibration_summary.csv\n")

## ---- Plots ----
df1 <- tibble(y = y26378, prob = p26378)
plot_score_dist(df1, "Locked GLM score distribution: GSE26378", file.path(dir_results, "F12_GSE26378_score_distribution.png"))
cur1 <- calibration_curve(y26378, p26378, bins = 10)
plot_calibration(cur1, "Calibration: GSE26378 (locked GLM)", file.path(dir_results, "F13_GSE26378_calibration.png"))

dfA <- tibble(y = yA, prob = p28750[keepA])
plot_score_dist(dfA, "Locked GLM score distribution: GSE28750(A)", file.path(dir_results, "F12_GSE28750A_score_distribution.png"))
curA <- calibration_curve(yA, p28750[keepA], bins = 6)
plot_calibration(curA, "Calibration: GSE28750(A) (locked GLM)", file.path(dir_results, "F13_GSE28750A_calibration.png"))

dfB <- tibble(y = yB, prob = p28750[keepB])
plot_score_dist(dfB, "Locked GLM score distribution: GSE28750(B)", file.path(dir_results, "F12_GSE28750B_score_distribution.png"))
curB <- calibration_curve(yB, p28750[keepB], bins = 8)
plot_calibration(curB, "Calibration: GSE28750(B) (locked GLM)", file.path(dir_results, "F13_GSE28750B_calibration.png"))

## ---- Exploratory thresholds (supplement only) ----
thr_ext_26378 <- find_thr_youden(y26378, p26378)
thr_ext_A <- find_thr_youden(yA, p28750[keepA])
thr_ext_B <- find_thr_youden(yB, p28750[keepB])

thr_sp90_26378 <- find_thr_spec(y26378, p26378, 0.90)
thr_sp90_A <- find_thr_spec(yA, p28750[keepA], 0.90)
thr_sp90_B <- find_thr_spec(yB, p28750[keepB], 0.90)

thr_tbl <- tibble(
  cohort = c("GSE26378", "GSE28750", "GSE28750"),
  endpoint = c(
    "Diagnosis (Septic_shock vs Control)",
    "Diagnosis A (Sepsis vs Post_surgical)",
    "Diagnosis B (Sepsis vs All controls)"
  ),
  threshold_train_locked = c(thr_train_youden, thr_train_youden, thr_train_youden),
  threshold_external_youden = c(thr_ext_26378, thr_ext_A, thr_ext_B),
  threshold_external_spec90 = c(thr_sp90_26378, thr_sp90_A, thr_sp90_B)
)

readr::write_csv(thr_tbl, file.path(dir_results, "T19_external_exploratory_thresholds.csv"))
cat("✅ Saved: 04_results/T19_external_exploratory_thresholds.csv\n")