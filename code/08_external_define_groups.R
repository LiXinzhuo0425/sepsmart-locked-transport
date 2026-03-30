## =========================================================
## 08_external_define_groups.R
## Purpose:
##   Define external cohort groups/outcomes (lock rules) and export flow counts
## Inputs:
##   02_data_clean/GSE26378_pheno.csv
##   02_data_clean/GSE28750_pheno.csv
## Outputs:
##   02_data_clean/GSE26378_pheno_labeled.csv
##   02_data_clean/GSE28750_pheno_labeled.csv
##   04_results/T11_external_flow_counts.csv
## =========================================================

source('03_scripts/00_paths.R')
options(stringsAsFactors = FALSE)

if (!requireNamespace('dplyr', quietly = TRUE)) install.packages('dplyr')
if (!requireNamespace('readr', quietly = TRUE)) install.packages('readr')
if (!requireNamespace('stringr', quietly = TRUE)) install.packages('stringr')

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
})

norm_txt <- function(x) {
  x <- as.character(x)
  x[is.na(x)] <- ''
  tolower(trimws(x))
}

# -------------------------
# 1) GSE26378 labeling
# -------------------------
p1 <- file.path(dir_clean, 'GSE26378_pheno.csv')
stopifnot(file.exists(p1))
ph1 <- readr::read_csv(p1, show_col_types = FALSE)

# group field (preferred): disease state:ch1
stopifnot('disease state:ch1' %in% colnames(ph1))
g1_raw <- norm_txt(ph1[['disease state:ch1']])

# Define diagnosis group: septic shock vs normal control
ph1 <- ph1 %>%
  mutate(
    cohort = 'GSE26378',
    dx_group = case_when(
      str_detect(g1_raw, 'septic') ~ 'Septic_shock',
      str_detect(g1_raw, 'control|normal|healthy') ~ 'Control',
      TRUE ~ NA_character_
    )
  )

# outcome field (preferred): outcome:ch1
if ('outcome:ch1' %in% colnames(ph1)) {
  o1_raw <- norm_txt(ph1[['outcome:ch1']])
  ph1 <- ph1 %>% mutate(
    mort = case_when(
      str_detect(o1_raw, 'non') & str_detect(o1_raw, 'surviv') ~ 1L,
      str_detect(o1_raw, 'surviv') ~ 0L,
      TRUE ~ NA_integer_
    )
  )
} else {
  ph1 <- ph1 %>% mutate(mort = NA_integer_)
}

# Apply minimal inclusion: must have dx_group; mort optional
ph1_labeled <- ph1 %>%
  mutate(include_dx = !is.na(dx_group))

readr::write_csv(ph1_labeled, file.path(dir_clean, 'GSE26378_pheno_labeled.csv'))

# -------------------------
# 2) GSE28750 labeling
# -------------------------
p2 <- file.path(dir_clean, 'GSE28750_pheno.csv')
stopifnot(file.exists(p2))
ph2 <- readr::read_csv(p2, show_col_types = FALSE)

stopifnot('health status:ch1' %in% colnames(ph2))
g2_raw <- norm_txt(ph2[['health status:ch1']])

ph2 <- ph2 %>%
  mutate(
    cohort = 'GSE28750',
    # original 3-class label
    dx_3class = case_when(
      str_detect(g2_raw, 'sepsis') ~ 'Sepsis',
      str_detect(g2_raw, 'healthy') ~ 'Healthy',
      str_detect(g2_raw, 'post') | str_detect(g2_raw, 'surg') ~ 'Post_surgical',
      TRUE ~ NA_character_
    ),
    # 2-class option A: Sepsis vs Post_surgical (inflammatory control)
    dx_A = case_when(
      dx_3class == 'Sepsis' ~ 'Sepsis',
      dx_3class == 'Post_surgical' ~ 'Control_PostSurg',
      TRUE ~ NA_character_
    ),
    # 2-class option B: Sepsis vs all non-sepsis (Healthy + Post_surgical)
    dx_B = case_when(
      dx_3class == 'Sepsis' ~ 'Sepsis',
      dx_3class %in% c('Healthy','Post_surgical') ~ 'Control_All',
      TRUE ~ NA_character_
    )
  )

ph2_labeled <- ph2 %>%
  mutate(
    include_A = !is.na(dx_A),
    include_B = !is.na(dx_B),
    mort = NA_integer_
  )

readr::write_csv(ph2_labeled, file.path(dir_clean, 'GSE28750_pheno_labeled.csv'))

# -------------------------
# 3) Flow counts table (T11)
# -------------------------
flow_26378 <- tibble(
  cohort='GSE26378',
  n_total=nrow(ph1_labeled),
  n_dx_included=sum(ph1_labeled$include_dx, na.rm=TRUE),
  n_dx_septic=sum(ph1_labeled$dx_group=='Septic_shock', na.rm=TRUE),
  n_dx_control=sum(ph1_labeled$dx_group=='Control', na.rm=TRUE),
  n_mort_available=sum(!is.na(ph1_labeled$mort), na.rm=TRUE),
  n_dead=sum(ph1_labeled$mort==1, na.rm=TRUE),
  n_alive=sum(ph1_labeled$mort==0, na.rm=TRUE)
)

flow_28750_A <- tibble(
  cohort='GSE28750 (A: Sepsis vs Post_surgical)',
  n_total=nrow(ph2_labeled),
  n_dx_included=sum(ph2_labeled$include_A, na.rm=TRUE),
  n_dx_septic=sum(ph2_labeled$dx_A=='Sepsis', na.rm=TRUE),
  n_dx_control=sum(ph2_labeled$dx_A=='Control_PostSurg', na.rm=TRUE),
  n_mort_available=0, n_dead=0, n_alive=0
)

flow_28750_B <- tibble(
  cohort='GSE28750 (B: Sepsis vs All controls)',
  n_total=nrow(ph2_labeled),
  n_dx_included=sum(ph2_labeled$include_B, na.rm=TRUE),
  n_dx_septic=sum(ph2_labeled$dx_B=='Sepsis', na.rm=TRUE),
  n_dx_control=sum(ph2_labeled$dx_B=='Control_All', na.rm=TRUE),
  n_mort_available=0, n_dead=0, n_alive=0
)

T11 <- bind_rows(flow_26378, flow_28750_A, flow_28750_B)
out_t11 <- file.path(dir_results, 'T11_external_flow_counts.csv')
readr::write_csv(T11, out_t11)

cat('✅ Done. Outputs:\n')
cat(' - ', file.path(dir_clean, 'GSE26378_pheno_labeled.csv'), '\n', sep='')
cat(' - ', file.path(dir_clean, 'GSE28750_pheno_labeled.csv'), '\n', sep='')
cat(' - ', out_t11, '\n', sep='')
