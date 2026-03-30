## =========================================================
## 01_check_env.R
## Environment check + install required packages (renv)
## =========================================================

source(file.path('03_scripts', '00_paths.R'))

options(renv.consent = TRUE)

## ---- 1) Basic R session info ----
msg('R version: ', R.version.string)
msg('Platform: ', R.version$platform)

## ---- 2) Ensure renv is active ----
if (!requireNamespace('renv', quietly = TRUE)) {
  install.packages('renv')
}
msg('renv version: ', as.character(utils::packageVersion('renv')))

## Activate project renv (safe if already active)
renv::activate(project_root)
msg('renv activated at: ', project_root)

## ---- 3) Install required packages via renv ----
cran_pkgs <- c(
  'data.table','dplyr','tibble','readr','stringr','ggplot2','pROC',
  'caret','glmnet','Matrix','patchwork','optparse','rmarkdown'
)

bioc_pkgs <- c(
  'BiocManager','GEOquery','limma','edgeR','DESeq2',
  'clusterProfiler','org.Hs.eg.db','AnnotationDbi','fgsea'
)

install_if_missing <- function(pkgs, repo = c('CRAN','Bioc')) {
  repo <- match.arg(repo)
  for (p in pkgs) {
    if (!requireNamespace(p, quietly = TRUE)) {
      msg('Installing ', p, ' (', repo, ') ...')
      if (repo == 'CRAN') {
        install.packages(p, dependencies = TRUE)
      } else {
        if (!requireNamespace('BiocManager', quietly = TRUE)) install.packages('BiocManager')
        BiocManager::install(p, update = FALSE, ask = FALSE)
      }
    } else {
      msg('OK: ', p, ' already installed')
    }
  }
}

install_if_missing(cran_pkgs, repo = 'CRAN')
install_if_missing(bioc_pkgs, repo = 'Bioc')

## ---- 4) Snapshot lockfile (records installed versions) ----
renv::snapshot(prompt = FALSE)
msg('renv snapshot complete. Lockfile updated: ', file.path(project_root, 'renv.lock'))

## ---- 5) Quick sanity checks ----
suppressPackageStartupMessages({
  library(GEOquery); library(limma); library(clusterProfiler); library(pROC)
})
msg('Loaded core packages successfully.')
msg('SepSMART genes: ', paste(genes_sepSMART, collapse = ', '))
