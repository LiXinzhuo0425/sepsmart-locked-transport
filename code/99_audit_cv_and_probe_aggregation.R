# ============================================================
# 03_scripts/99_audit_cv_and_probe_aggregation.R
# Audit: CV folds / repeats / probe aggregation rules
# Author: (auto)
# ============================================================

# ---------- 0) safety ----------
options(stringsAsFactors = FALSE)

project_root <- getwd()
cat("[INFO] project_root =", project_root, "\n")

scripts_dir <- file.path(project_root, "03_scripts")
if (!dir.exists(scripts_dir)) stop("Not found: ", scripts_dir)

# try load your path config if exists
paths_r <- file.path(scripts_dir, "00_paths.R")
if (file.exists(paths_r)) {
  try(source(paths_r), silent = TRUE)
  cat("[INFO] Sourced 00_paths.R\n")
}

# ---------- 1) helpers ----------
read_lines_safe <- function(f) {
  x <- try(readLines(f, warn = FALSE, encoding = "UTF-8"), silent = TRUE)
  if (inherits(x, "try-error")) {
    x <- try(readLines(f, warn = FALSE), silent = TRUE)
  }
  if (inherits(x, "try-error")) return(character(0))
  x
}

grep_hits <- function(lines, patterns, ignore_case = TRUE) {
  # returns data.frame with line_no + line_text + matched_pattern
  out <- list()
  for (pat in patterns) {
    idx <- grep(pat, lines, ignore.case = ignore_case, perl = TRUE)
    if (length(idx)) {
      out[[length(out) + 1]] <- data.frame(
        line_no = idx,
        line_text = lines[idx],
        pattern = rep(pat, length(idx)),
        stringsAsFactors = FALSE
      )
    }
  }
  if (!length(out)) return(data.frame())
  do.call(rbind, out)
}

pretty_print_hits <- function(df, file, max_lines = 200) {
  if (!nrow(df)) return(invisible(NULL))
  df <- df[order(df$line_no), , drop = FALSE]
  if (nrow(df) > max_lines) df <- df[1:max_lines, , drop = FALSE]
  cat("\n------------------------------------------------------------\n")
  cat("FILE:", file, "\n")
  cat("HITS:", nrow(df), "\n")
  cat("------------------------------------------------------------\n")
  for (i in seq_len(nrow(df))) {
    cat(sprintf("L%-5d  %s\n", df$line_no[i], df$line_text[i]))
  }
  invisible(NULL)
}

# ---------- 2) scan all R scripts ----------
r_files <- list.files(scripts_dir, pattern = "\\.R$", full.names = TRUE, recursive = TRUE)
if (!length(r_files)) stop("No .R files under: ", scripts_dir)

cat("[INFO] Found", length(r_files), "R scripts under", scripts_dir, "\n")

# ---------- 3) patterns to detect CV config ----------
# We intentionally include many frameworks:
# - caret: trainControl(method="cv"/"repeatedcv", number=, repeats=)
# - glmnet: cv.glmnet(..., nfolds=)
# - rsample: vfold_cv(v=), bootstraps(times=)
# - custom folds: createFolds, groupKFold
cv_patterns <- c(
  # caret
  "trainControl\\s*\\(",
  "method\\s*=\\s*['\"](cv|repeatedcv|LOOCV|boot|boots632)['\"]",
  "number\\s*=\\s*\\d+",
  "repeats\\s*=\\s*\\d+",
  # glmnet
  "cv\\.glmnet\\s*\\(",
  "nfolds\\s*=\\s*\\d+",
  # rsample / tidymodels
  "vfold_cv\\s*\\(",
  "\\bv\\s*=\\s*\\d+",
  "bootstraps\\s*\\(",
  "times\\s*=\\s*\\d+",
  # custom folds
  "createFolds\\s*\\(",
  "groupKFold\\s*\\(",
  "KFold|kfold|k-fold",
  "\\bfold\\b",
  "set\\.seed\\s*\\("
)

# ---------- 4) patterns to detect probe aggregation / mapping rules ----------
# This captures:
# - limma::avereps
# - WGCNA::collapseRows
# - aggregate/summarise
# - choose max/mean/median
# - “probe”, “Gene Symbol”, “fData” mapping
probe_patterns <- c(
  "avereps\\s*\\(",
  "collapseRows\\s*\\(",
  "aggregate\\s*\\(",
  "rowsum\\s*\\(",
  "summaris(e|z)\\s*\\(",
  "group_by\\s*\\(",
  "distinct\\s*\\(",
  "duplicated\\s*\\(",
  "probe",
  "probeset|probeset_id|probe_id|ID_REF",
  "Gene\\s*Symbol|gene_symbol|symbol",
  "fData\\s*\\(",
  "annotation|annot",
  "mean\\s*\\(|median\\s*\\(|max\\s*\\(",
  "which\\.max\\s*\\(",
  "select\\s*\\(|slice_max\\s*\\(",
  "largest|highest|variance|IQR|sd\\s*\\("
)

# ---------- 5) run scan ----------
cv_report <- list()
probe_report <- list()

for (f in r_files) {
  lines <- read_lines_safe(f)
  if (!length(lines)) next
  
  cv_hits <- grep_hits(lines, cv_patterns)
  if (nrow(cv_hits)) {
    cv_hits$file <- f
    cv_report[[length(cv_report) + 1]] <- cv_hits
  }
  
  pr_hits <- grep_hits(lines, probe_patterns)
  if (nrow(pr_hits)) {
    pr_hits$file <- f
    probe_report[[length(probe_report) + 1]] <- pr_hits
  }
}

cv_df <- if (length(cv_report)) do.call(rbind, cv_report) else data.frame()
pr_df <- if (length(probe_report)) do.call(rbind, probe_report) else data.frame()

# ---------- 6) summarize likely CV settings ----------
infer_cv_summary <- function(cv_df) {
  if (!nrow(cv_df)) {
    return(list(
      caret_method = NA, caret_number = NA, caret_repeats = NA,
      glmnet_nfolds = NA
    ))
  }
  txt <- paste(cv_df$line_text, collapse = "\n")
  
  # caret method
  m <- regmatches(txt, regexpr("method\\s*=\\s*['\"][^'\"]+['\"]", txt, perl = TRUE))
  caret_method <- if (length(m) && nchar(m)) gsub(".*['\"]|['\"].*", "", m) else NA
  
  # caret number
  n <- regmatches(txt, regexpr("number\\s*=\\s*\\d+", txt, perl = TRUE))
  caret_number <- if (length(n) && nchar(n)) as.integer(gsub(".*=\\s*", "", n)) else NA
  
  # caret repeats
  r <- regmatches(txt, regexpr("repeats\\s*=\\s*\\d+", txt, perl = TRUE))
  caret_repeats <- if (length(r) && nchar(r)) as.integer(gsub(".*=\\s*", "", r)) else NA
  
  # glmnet nfolds
  g <- regmatches(txt, regexpr("nfolds\\s*=\\s*\\d+", txt, perl = TRUE))
  glmnet_nfolds <- if (length(g) && nchar(g)) as.integer(gsub(".*=\\s*", "", g)) else NA
  
  list(
    caret_method = caret_method,
    caret_number = caret_number,
    caret_repeats = caret_repeats,
    glmnet_nfolds = glmnet_nfolds
  )
}

cv_summary <- infer_cv_summary(cv_df)

# ---------- 7) write audit outputs ----------
out_dir <- if (exists("results_dir")) {
  get("results_dir")
} else {
  file.path(project_root, "04_results")
}
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

out_csv_cv  <- file.path(out_dir, "T99_audit_CV_hits.csv")
out_csv_pr  <- file.path(out_dir, "T99_audit_probe_hits.csv")
out_txt_sum <- file.path(out_dir, "T99_audit_summary.txt")

if (nrow(cv_df)) write.csv(cv_df, out_csv_cv, row.names = FALSE)
if (nrow(pr_df)) write.csv(pr_df, out_csv_pr, row.names = FALSE)

sum_lines <- c(
  "=== Audit summary: CV folds / repeats / probe aggregation ===",
  paste0("project_root: ", project_root),
  paste0("scripts_dir : ", scripts_dir),
  "",
  "---- CV summary (best-effort inference from hits) ----",
  paste0("caret method   : ", cv_summary$caret_method),
  paste0("caret number   : ", cv_summary$caret_number),
  paste0("caret repeats  : ", cv_summary$caret_repeats),
  paste0("glmnet nfolds  : ", cv_summary$glmnet_nfolds),
  "",
  "---- Output files ----",
  paste0("CV hits   : ", out_csv_cv),
  paste0("Probe hits: ", out_csv_pr),
  "",
  "Notes:",
  "- If caret_method == 'repeatedcv', then repeats indicates repeated K-fold CV.",
  "- If glmnet_nfolds is present, that indicates K-fold CV inside cv.glmnet.",
  "- Probe aggregation rules must be read from the specific lines in T99_audit_probe_hits.csv"
)
writeLines(sum_lines, out_txt_sum)

# ---------- 8) print to console (human reading) ----------
cat("\n================= AUDIT SUMMARY =================\n")
cat(paste(sum_lines, collapse = "\n"), "\n")
cat("=================================================\n")

# print detailed hits (top)
if (nrow(cv_df)) {
  cat("\n\n### CV-related hits (grouped by file)\n")
  for (ff in unique(cv_df$file)) {
    pretty_print_hits(cv_df[cv_df$file == ff, c("line_no","line_text","pattern")], ff, max_lines = 120)
  }
} else {
  cat("\n[WARN] No CV-related hits found in scripts.\n")
}

if (nrow(pr_df)) {
  cat("\n\n### Probe/mapping-related hits (grouped by file)\n")
  for (ff in unique(pr_df$file)) {
    pretty_print_hits(pr_df[pr_df$file == ff, c("line_no","line_text","pattern")], ff, max_lines = 160)
  }
} else {
  cat("\n[WARN] No probe/mapping-related hits found in scripts.\n")
}

cat("\n✅ Done. Review these files:\n- ", out_txt_sum,
    "\n- ", out_csv_cv,
    "\n- ", out_csv_pr, "\n", sep = "")