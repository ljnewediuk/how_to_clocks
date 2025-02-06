
# EWAS to find CpGs correlated with sample characteristics. Remove any probes 
# correlated (p < 0.05) with tissue, sex, and population andkeep top ~ 10% of 
# probes correlated with age. 

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
source('functions/groupLimma.R')

# 2 - Load data ===

# Normalized betas and sample info
meth_dat <- readRDS('input/norm_betas.rds') 

# 3 - Run EWAS ====

# Fit limmas for variable of interest, get CpGs correlated with age based on 
# adjusted p-value, and extract all variables correlated with age across the 
# variable

# TISSUE 
# (results in a very small number of CpG sites common to all tissues)
tissue_list <- unique(meth_dat$Spec)
for(i in 1:length(tissue_list)) {
  # Fit limma
  lim_cgs <- groupLimma(dat = meth_dat, variable = 'Spec', value = tissue_list[i],
                        control_vars <- c('sex', 'Population'), pval = 1e-3)
  # Get CpG sites correlated with age in all categories
  if(i == 1) tissue_cgs <- lim_cgs
  if(i > 1) tissue_cgs <- tissue_cgs[tissue_cgs %in% lim_cgs]
}

# SEX 
# (results in an appropriate number of CpG sites)
sex_list <- c('M', 'F')
for(i in 1:length(sex_list)) {
  # Fit limma (can't control for population because one population is all male)
  lim_cgs <- groupLimma(dat = meth_dat, variable = 'sex', value = sex_list[i],
                        control_vars <- c('Spec'), pval = 1e-3)
  # Get CpG sites correlated with age in all categories
  if(i == 1) sex_cgs <- lim_cgs
  if(i > 1) sex_cgs <- sex_cgs[sex_cgs %in% lim_cgs]
}

# POPULATION 
# (no CpG sites are common across all populations; suggests we might need to be 
# more conservative in excluding sites)
pop_list <- unique(meth_dat$Population)
for(i in 1:length(pop_list)) {
  # Can't control for tissue by population because most populations are all
  # the same tissue; can't control for sex in NW population because all male
  if(! pop_list[i] == 'NW') {
    control_vars <- 'sex'
  } else {
    control_vars <- NULL
  }
  # Get CpG sites correlated with age in all categories
  lim_cgs <- groupLimma(dat = meth_dat, variable = 'Population', value = pop_list[i],
                        control_vars = control_vars, pval = 1e-3)
  if(i == 1) pop_cgs <- lim_cgs
  if(i > 1) pop_cgs <- pop_cgs[pop_cgs %in% lim_cgs]
}

# TISSUE & SEX
# Sites related to age for all tissues and sexes
tissue_sex_cgs <- tissue_cgs[tissue_cgs %in% sex_cgs]

# 6 - Save ====

saveRDS(tissue_cgs, 'output/meth_ewas_tis.rds')
saveRDS(sex_cgs, 'output/meth_ewas_sex.rds')
saveRDS(tissue_sex_cgs, 'output/meth_ewas_tis_sex.rds')
