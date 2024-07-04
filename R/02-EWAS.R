
# 02 - EWAS to find CpGs correlated with sample characteristics

# Remove any probes correlated (p < 0.05) with tissue, sex, and population and 
# keep top ~ 10% of probes correlated with age. 

library(tidyverse)

#   1) Two EWAS on age with blood and skin separated
#   2) 11 EWAS on age with all populations separated
#   3) Two EWAS on age with males and females separated
#   4) One EWAS on tissue, age, sex, and population, adjusting p-value for 
#      either sex, tissue, or population

# The design matrix needs to be formatted with samples in rows and coefficients 
# in columns

# 1 - Prep workspace ====

# Load EWAS function
source('functions/fitLimma.R')

# 2 - Load data ===

# Normalized betas and sample info
meth_dat <- readRDS('input/norm_betas.rds') 

# 3 - Subset data by tissue, population, and sex ====

# Split by population
split_pop <- meth_dat %>%
  group_by(Population) %>%
  group_split()

# Split by tissue
split_spec <- meth_dat %>%
  group_by(Spec) %>%
  group_split()

# Split by sex
split_sex <- meth_dat %>%
  group_by(sex) %>%
  group_split()

# 4 - Run EWAS ====

full_limma <- fitLimma(meth_dat = meth_dat)

# Add adjusted pvalues for variable of interest

pvals_full <- full_limma %>%
  mutate(across(matches('pval'), list(adj = function(x) p.adjust(x, method = 'BH'))))

# 5 - Correlated positions ====

# Get pvlaues either > 0.05 or < 0.05 (depending on whether we want the Cgs
# correlated or uncorrelated with the variable)

# Tissue (not correlated with)
cgs_cor_w_tissue <- pvals_full %>%
  select(CGid, matches('pval_adj') & matches('Spec') & ! 
           matches('Intercept') & ! matches('age')) %>%
  filter(if_any(everything(), ~ . < 0.05)) %>%
  pull(CGid)

# Sex (not correlated with)
cgs_cor_w_sex <- pvals_full %>%
  select(CGid, matches('pval_adj') & matches('sex') & ! 
           matches('Intercept') & ! matches('age')) %>%
  filter(if_any(everything(), ~ . < 0.05)) %>%
  pull(CGid)

# Population (not correlated with)
cgs_cor_w_pop <- pvals_full %>%
  select(CGid, matches('pval_adj') & matches('Population') & ! 
           matches('Intercept') & ! matches('age')) %>%
  filter(if_any(everything(), ~ . < 0.05)) %>%
  pull(CGid)

# Age (correlated with)
cgs_cor_w_age <- pvals_full %>%
  select(CGid, matches('pval_adj') & matches('age') & ! matches('Intercept')) %>%
  filter(if_any(everything(), ~ . < 0.05)) %>%
  pull(CGid)

# 6 - Save new data ====

# Reduced Cg positions based on lm
meth_dat_red <- meth_dat %>% 
  select(cgs_cor_w_age) %>%
  select(! any_of(c(cgs_cor_w_pop, cgs_cor_w_sex, cgs_cor_w_tissue)))

# Add back sample info
new_meth_dat <- meth_dat %>%
  select(sampleId:Population) %>%
  bind_cols(meth_dat_red)

saveRDS(new_meth_dat, 'output/lm_norm_betas.rds')














# 5 Run limma for EWAS on age separately for blood and skin ====

# Blood/age
blood_age_limma <- fit_limma(betas = sample_sheet_blood, 
                             mod_formula = c('age', 'sex'),
                             relatives = relatives_list[[i]]) %>%
  mutate(age_pval_adj = p.adjust(age_pval, method = 'BH')) %>%
  mutate(sig = ifelse(age_pval < 10e-6, 'yes', 'no'))

# Skin/age
skin_age_limma <- fit_limma(betas = sample_sheet_skin, 
                            mod_formula = c('age', 'sex'),
                            relatives = relatives_list[[i]]) %>%
  mutate(age_pval_adj = p.adjust(age_pval, method = 'BH')) %>%
  mutate(sig = ifelse(age_pval < 10e-6, 'yes', 'no'))

# Make table of p values correlated with age
print(
  data.frame(
    tissue = c(rep('skin', times = 3), rep('blood', times = 3)),
    pval = rep(c('10e-6', '10e-7', '10e-8'), times = 2),
    n = c(nrow(skin_age_limma[skin_age_limma$age_pval < 10e-6,]),
          nrow(skin_age_limma[skin_age_limma$age_pval < 10e-7,]),
          nrow(skin_age_limma[skin_age_limma$age_pval < 10e-8,]),
          nrow(blood_age_limma[blood_age_limma$age_pval < 10e-6,]),
          nrow(blood_age_limma[blood_age_limma$age_pval < 10e-7,]),
          nrow(blood_age_limma[blood_age_limma$age_pval < 10e-8,])))
)

# Use all Cgs highly correlated with age (p < 10e-6) and exclude Cgs correlated 
# (p < 0.05) with sex
Cgs_sample <- skin_age_limma %>%
  # Cgs correlated with age
  rbind(blood_age_limma) %>%
  filter(age_pval_adj < 10e-6) %>%
  # Remove Cgs correlated with sex and tissue
  filter(! CGid %in% c(cgs_cor_w_sex)) %>%
  pull(CGid)  %>% 
  unique()

# Assign Cgs as either with blocking for relatives or without
assign(paste0('Cgs_sample_', names(relatives_list)[i]), Cgs_sample)



# 6 Difference between Cgs selected with vs. without blocking by relatives ====

setdiff(Cgs_sample_no_relatives, Cgs_sample_relatives)

# Result:
# 8 Cpgs different: "cg20167048" "cg05631094" "cg06074849" "cg08383062" 
# "cg15148667" "cg16787065" "cg22661206" "cg08965235"
# 
# This is a ~0.2% difference from the feature set we used in the original clock,
# and ~0.02% of all Cpgs in the array.

# 7 Save ====
