## =========================================================
## 11_train_fixed_dx_model_on_GSE65682.R (v3 FIXED)
## Fix v3:
##   - Harmonize gene symbols (aliases -> canonical), e.g., C19ORF59 -> MCEMP1
##   - Allow training with >=5 genes if a gene is truly absent in platform
##   - Export expression with all 6 genes (missing filled with NA) for pipeline stability
## =========================================================

source('03_scripts/00_paths.R')
options(stringsAsFactors = FALSE)

SEPSMART_GENES <- c('RETN','MCEMP1','CYP1B1','S100A12','S100A8','HK3')
MIN_GENES_FOR_TRAIN <- 5

pkgs <- c('GEOquery','Biobase','dplyr','readr','stringr','tibble','data.table','glmnet','pROC','ggplot2')
missing <- pkgs[!vapply(pkgs, requireNamespace, quietly=TRUE, FUN.VALUE=logical(1))]
if (length(missing)>0) {
  if (requireNamespace('renv', quietly=TRUE)) try(renv::install(missing), silent=TRUE)
  missing2 <- missing[!vapply(missing, requireNamespace, quietly=TRUE, FUN.VALUE=logical(1))]
  if (length(missing2)>0) install.packages(missing2)
}
suppressPackageStartupMessages({
  library(GEOquery); library(Biobase)
  library(dplyr); library(readr); library(stringr); library(tibble)
  library(data.table); library(glmnet); library(pROC); library(ggplot2)
})

timestamp <- function() format(Sys.time(), '%Y-%m-%d %H:%M:%S')
safe_mkdir <- function(d) if (!dir.exists(d)) dir.create(d, recursive=TRUE, showWarnings=FALSE)

norm_symbol_vec <- function(x) {
  x <- as.character(x); x[is.na(x)] <- ''
  x <- toupper(x)
  x <- str_replace_all(x, '///', ';')
  x <- str_replace_all(x, '\\|', ';')
  x <- str_replace_all(x, ',', ';')
  x <- str_replace_all(x, '\\s+', '')
  x
}

unwrap_eset <- function(gobj) {
  if (inherits(gobj,'ExpressionSet')) return(gobj)
  if (is.list(gobj) && length(gobj)>=1) {
    for (i in seq_along(gobj)) if (inherits(gobj[[i]],'ExpressionSet')) return(gobj[[i]])
  }
  NULL
}

# --- Alias map: extend if needed ---
ALIAS_MAP <- c(
  'C19ORF59'='MCEMP1',
  'S100A8A'='S100A8'  # rare, placeholder
)

harmonize_gene_symbols <- function(gene_mat) {
  rn <- toupper(rownames(gene_mat))
  rn2 <- rn
  for (a in names(ALIAS_MAP)) rn2[rn2 == a] <- ALIAS_MAP[[a]]
  rownames(gene_mat) <- rn2
  # if duplicates created by aliasing, collapse by median
  if (any(duplicated(rn2))) {
    dt <- data.table::as.data.table(gene_mat, keep.rownames='symbol')
    cols <- setdiff(colnames(dt), 'symbol')
    dt2 <- dt[, lapply(.SD, median, na.rm=TRUE), by=symbol, .SDcols=cols]
    m <- as.matrix(dt2[, ..cols])
    rownames(m) <- dt2$symbol
    return(m)
  }
  gene_mat
}

probe_symbol_map <- function(eset) {
  fd <- tryCatch(Biobase::fData(eset), error=function(e) NULL)
  if (!is.null(fd) && nrow(fd)>0) {
    cn <- colnames(fd)
    sym_cols <- cn[str_detect(tolower(cn), 'gene.?symbol|symbol|genesymbol')]
    if (length(sym_cols)>0) {
      scores <- sapply(sym_cols, function(cc) sum(!is.na(fd[[cc]]) & fd[[cc]]!=''))
      sym_col <- sym_cols[which.max(scores)]
      m <- data.frame(probe=rownames(fd), symbol=norm_symbol_vec(fd[[sym_col]]), stringsAsFactors=FALSE)
      return(list(map=m, note=paste0('fData:', sym_col)))
    }
  }

  gpl_id <- tryCatch(Biobase::annotation(eset), error=function(e) NA_character_)
  if (is.na(gpl_id) || !nzchar(gpl_id)) return(list(map=data.frame(), note='No fData symbol; no GPL id'))
  gpl_obj <- tryCatch(GEOquery::getGEO(gpl_id), error=function(e) NULL)
  if (is.null(gpl_obj)) return(list(map=data.frame(), note=paste0('Fail GPL download: ', gpl_id)))
  tab <- tryCatch(GEOquery::Table(gpl_obj), error=function(e) NULL)
  if (is.null(tab) || nrow(tab)==0) return(list(map=data.frame(), note=paste0('Empty GPL table: ', gpl_id)))

  id_col <- if ('ID' %in% colnames(tab)) 'ID' else colnames(tab)[1]
  cn <- colnames(tab)
  sym_cols <- cn[str_detect(tolower(cn), 'gene.?symbol|symbol|gene symbol')]
  if (length(sym_cols)==0) return(list(map=data.frame(), note=paste0('No symbol col in GPL: ', gpl_id)))
  scores <- sapply(sym_cols, function(cc) sum(!is.na(tab[[cc]]) & tab[[cc]]!=''))
  sym_col <- sym_cols[which.max(scores)]

  m <- tab[, c(id_col, sym_col)]
  colnames(m) <- c('probe','symbol')
  m$symbol <- norm_symbol_vec(m$symbol)
  m <- as.data.frame(m, stringsAsFactors=FALSE)
  list(map=m, note=paste0('GPL:', gpl_id, '/', sym_col))
}

collapse_to_gene <- function(expr_mat, map_df) {
  if (is.null(map_df) || nrow(map_df)==0) stop('Empty probe->symbol map; cannot collapse.')
  map_df <- as.data.frame(map_df, stringsAsFactors=FALSE)
  map_df$probe <- as.character(map_df$probe)
  map_df$symbol <- norm_symbol_vec(map_df$symbol)
  map_df <- map_df[!is.na(map_df$probe) & map_df$probe!='' & !is.na(map_df$symbol) & map_df$symbol!='', , drop=FALSE]

  sym_split <- str_split(map_df$symbol, ';')
  map_exp <- data.frame(
    probe = rep(map_df$probe, lengths(sym_split)),
    symbol = unlist(sym_split, use.names=FALSE),
    stringsAsFactors = FALSE
  )
  map_exp$symbol <- toupper(map_exp$symbol)
  map_exp <- map_exp[!is.na(map_exp$symbol) & map_exp$symbol!='', , drop=FALSE]
  map_exp <- unique(map_exp)
  map_exp <- map_exp[map_exp$probe %in% rownames(expr_mat), , drop=FALSE]
  if (nrow(map_exp)==0) stop('No overlapping probes between expression and mapping.')

  dt <- data.table::as.data.table(expr_mat, keep.rownames='probe')
  dt <- merge(data.table::as.data.table(map_exp), dt, by='probe', all=FALSE)
  cols <- setdiff(colnames(dt), c('probe','symbol'))
  gene_dt <- dt[, lapply(.SD, median, na.rm=TRUE), by=symbol, .SDcols=cols]
  gene_mat <- as.matrix(gene_dt[, ..cols])
  rownames(gene_mat) <- gene_dt$symbol
  gene_mat
}

# -------------------------
# Load pheno
# -------------------------
ph_path <- file.path(dir_clean, 'GSE65682_single_platform_pheno.csv')
stopifnot(file.exists(ph_path))
ph <- readr::read_csv(ph_path, show_col_types=FALSE) %>% mutate(.sid=as.character(geo_accession))
stopifnot('abdominal_sepsis_and_controls:ch1' %in% names(ph))

# -------------------------
# Load / Download ExpressionSet
# -------------------------
raw_dir <- file.path(dir_raw, 'GSE65682')
safe_mkdir(raw_dir)
eset_path <- file.path(raw_dir, 'GSE65682_eset.rds')

if (file.exists(eset_path)) {
  message('[', timestamp(), '] Loading cached eset: ', eset_path)
  eset <- readRDS(eset_path)
} else {
  message('[', timestamp(), '] Downloading GSE65682 via GEOquery::getGEO ...')
  gobj <- GEOquery::getGEO('GSE65682', GSEMatrix=TRUE, getGPL=TRUE)
  eset <- unwrap_eset(gobj)
  if (is.null(eset)) stop('Cannot unwrap getGEO output to ExpressionSet for GSE65682')
  saveRDS(eset, eset_path)
  message('[', timestamp(), '] Saved eset: ', eset_path)
}

expr <- Biobase::exprs(eset)  # probe x sample

# -------------------------
# Map -> gene, harmonize symbols, then extract SepSMART genes
# -------------------------
mp <- probe_symbol_map(eset)
gene_mat <- collapse_to_gene(expr, mp$map)  # gene x sample
gene_mat <- harmonize_gene_symbols(gene_mat)

genes_have <- intersect(SEPSMART_GENES, rownames(gene_mat))
genes_miss <- setdiff(SEPSMART_GENES, genes_have)

if (length(genes_have) < MIN_GENES_FOR_TRAIN) {
  stop('GSE65682: too many SepSMART genes missing. have=', length(genes_have), ' missing=', paste(genes_miss, collapse=';'))
}

if (length(genes_miss) > 0) {
  warning('GSE65682: SepSMART gene(s) missing after alias harmonization: ', paste(genes_miss, collapse=';'),
          '. Proceeding with available genes (>=', MIN_GENES_FOR_TRAIN, ').')
}

# Build 6-gene matrix with NA rows for missing genes (for stable downstream I/O)
emat6 <- matrix(NA_real_, nrow=length(SEPSMART_GENES), ncol=ncol(gene_mat))
rownames(emat6) <- SEPSMART_GENES
colnames(emat6) <- colnames(gene_mat)
emat6[genes_have, ] <- gene_mat[genes_have, , drop=FALSE]

# Export WITH Gene column
out_expr_path <- file.path(dir_clean, 'GSE65682_SepSMART_expr.csv')
out_df <- data.frame(Gene=SEPSMART_GENES, emat6, check.names=FALSE)
readr::write_csv(out_df, out_expr_path)
message('[', timestamp(), '] Saved 6-gene expression: ', out_expr_path, ' (mapping=', mp$note, ')')

# -------------------------
# Align with pheno
# -------------------------
common <- intersect(colnames(emat6), ph$.sid)
if (length(common) < 30) stop('Too few overlapping samples between expression and pheno. common=', length(common))
ph2 <- ph %>% filter(.sid %in% common) %>% arrange(match(.sid, common))
emat6 <- emat6[, common, drop=FALSE]

# -------------------------
# Define training subset: abdo_s vs ctrl_GI
# -------------------------
grp <- as.character(ph2[['abdominal_sepsis_and_controls:ch1']])
keep <- !is.na(grp) & grp %in% c('abdo_s','ctrl_GI')
ph_tr <- ph2[keep, ]
emat_tr <- emat6[, keep, drop=FALSE]
y <- ifelse(ph_tr[['abdominal_sepsis_and_controls:ch1']]=='abdo_s', 1L, 0L)
cat('[INFO] Training subset n=', length(y), ' (abdo_s=', sum(y==1), ', ctrl_GI=', sum(y==0), ')\n', sep='')

# Use only available genes for training (drop missing)
train_genes <- intersect(SEPSMART_GENES, rownames(emat_tr))
train_genes <- train_genes[!is.na(rowMeans(emat_tr[train_genes, , drop=FALSE]))]  # remove all-NA rows
if (length(train_genes) < MIN_GENES_FOR_TRAIN) stop('After NA removal, too few genes for training: ', length(train_genes))
cat('[INFO] Genes used for training: ', paste(train_genes, collapse=','), '\n', sep='')

# -------------------------
# Standardize using TRAINING data and lock mean/sd
# -------------------------
X <- t(emat_tr[train_genes, , drop=FALSE])  # samples x genes
mu <- colMeans(X, na.rm=TRUE)
sdv <- apply(X, 2, sd, na.rm=TRUE)
sdv[sdv==0 | is.na(sdv)] <- 1
Xz <- scale(X, center=mu, scale=sdv)

# GLM
df <- data.frame(y=y, Xz)
glm_fit <- glm(y ~ ., data=df, family=binomial())

# LASSO
set.seed(1)
cvfit <- cv.glmnet(as.matrix(Xz), y, family='binomial', alpha=1, nfolds=5, type.measure='auc')
beta_min <- as.matrix(coef(cvfit, s='lambda.min'))
beta_1se <- as.matrix(coef(cvfit, s='lambda.1se'))

# Internal 5-fold CV predictions for GLM
set.seed(1)
k <- 5
foldid <- sample(rep(1:k, length.out=nrow(Xz)))
pred_glm <- rep(NA_real_, length(y))
for (f in 1:k) {
  tr <- which(foldid != f); te <- which(foldid == f)
  dtr <- data.frame(y=y[tr], Xz[tr, , drop=FALSE])
  fit <- glm(y ~ ., data=dtr, family=binomial())
  pred_glm[te] <- predict(fit, newdata=data.frame(Xz[te, , drop=FALSE]), type='response')
}

roc_glm <- pROC::roc(y, pred_glm, quiet=TRUE)
auc_glm <- as.numeric(pROC::auc(roc_glm))
ci_glm <- as.numeric(pROC::ci.auc(roc_glm))

dfroc <- data.frame(fpr=1-roc_glm$specificities, tpr=roc_glm$sensitivities)
p <- ggplot(dfroc, aes(x=fpr, y=tpr)) +
  geom_line(linewidth=1) + geom_abline(slope=1, intercept=0, linetype=2) +
  theme_bw() + labs(title=sprintf('GSE65682 internal CV ROC (GLM) AUC=%.3f', auc_glm), x='False Positive Rate', y='True Positive Rate')
ggsave(file.path(dir_results,'F10_GSE65682_dx_internal_ROC.png'), p, width=6.5, height=5, dpi=220)

# Export scaling params (include only genes actually used)
T15 <- tibble(gene=train_genes, mean=as.numeric(mu[train_genes]), sd=as.numeric(sdv[train_genes]))
readr::write_csv(T15, file.path(dir_results,'T15_GSE65682_dx_scaling_params.csv'))

# Export coefficients
coef_glm <- coef(glm_fit)
coef_tbl_glm <- tibble(model='GLM', term=names(coef_glm), coefficient=as.numeric(coef_glm))
coef_tbl_min <- tibble(model='LASSO_min', term=rownames(beta_min), coefficient=as.numeric(beta_min[,1]))
coef_tbl_1se <- tibble(model='LASSO_1se', term=rownames(beta_1se), coefficient=as.numeric(beta_1se[,1]))
T14 <- bind_rows(coef_tbl_glm, coef_tbl_min, coef_tbl_1se)
readr::write_csv(T14, file.path(dir_results,'T14_GSE65682_dx_model_coefficients.csv'))

# Export internal CV predictions
T14b <- tibble(sample_id=ph_tr$.sid, y=y, pred_glm_cv=pred_glm)
readr::write_csv(T14b, file.path(dir_results,'T14b_GSE65682_dx_internal_cv_predictions.csv'))

cat('✅ Saved: 02_data_clean/GSE65682_SepSMART_expr.csv\n')
cat('✅ Saved: T14 coefficients, T15 scaling, F10 ROC, T14b predictions.\n')
cat('   Internal CV AUC(GLM)=', sprintf('%.3f', auc_glm), ' (95%CI ', sprintf('%.3f', ci_glm[1]), '-', sprintf('%.3f', ci_glm[3]), ')\n', sep='')
cat('   Mapping used: ', mp$note, '\n', sep='')
if (length(genes_miss)>0) cat('   Missing genes (reported): ', paste(genes_miss, collapse=';'), '\n', sep='')
