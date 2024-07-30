
# 02 - EWAS to find CpGs correlated with sample characteristics

# Remove any probes correlated (p < 0.05) with tissue, sex, and population and 
# keep top ~ 10% of probes correlated with age. 

#   1) Two EWAS on age with blood and skin separated
#   2) 11 EWAS on age with all populations separated
#   3) Two EWAS on age with males and females separated
#   4) One EWAS on tissue, age, sex, and population, adjusting p-value for 
#      either sex, tissue, or population

# The design matrix needs to be formatted with samples in rows and coefficients 
# in columns

# 1 - Prep workspace ====

# Load libraries
library(tidyverse)

# Load EWAS function
source('functions/fitLimma.R')

# 2 - Load data ===

# Normalized betas and sample info
meth_dat <- readRDS('input/norm_betas.rds') 

# 3 - Run EWAS ====

# Fit lm to test probe relationships with age, sex, tissue, and population
# while also controlling for the other variables
full_limma <- fitLimma(meth_dat = meth_dat)

# Add adjusted pvalues for variable of interest

pvals_full <- full_limma %>%
  mutate(across(matches('pval'), list(adj = function(x) p.adjust(x, method = 'BH'))))

# 4 - Correlated positions ====

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

# 5 - Select positions ====

# Reduced Cg positions based on lm (remove all correlated with sex, tissue, and
# population and keep remaining Cgs correlated with age)
meth_dat_full <- meth_dat %>% 
  select(c(sampleId:Population, cgs_cor_w_age)) %>%
  select(! any_of(c(cgs_cor_w_pop, cgs_cor_w_sex, cgs_cor_w_tissue)))

# Cg positions not correlated with population, but correlated with age
meth_dat_pop <- meth_dat %>% 
  select(c(sampleId:Population, cgs_cor_w_age)) %>%
  select(! any_of(cgs_cor_w_pop))

# Cg positions not correlated with tissue, but correlated with age
meth_dat_tissue <- meth_dat %>% 
  select(c(sampleId:Population, cgs_cor_w_age)) %>%
  select(! any_of(cgs_cor_w_tissue))

# Cg positions not correlated with sex, but correlated with age
meth_dat_sex <- meth_dat %>% 
  select(c(sampleId:Population, cgs_cor_w_age)) %>%
  select(! any_of(cgs_cor_w_sex))

# 6 - Save ====

saveRDS(meth_dat_full, 'output/meth_ewas_all.rds')
saveRDS(meth_dat_pop, 'output/meth_ewas_pop.rds')
saveRDS(meth_dat_tissue, 'output/meth_ewas_tis.rds')
saveRDS(meth_dat_sex, 'output/meth_ewas_sex.rds')
