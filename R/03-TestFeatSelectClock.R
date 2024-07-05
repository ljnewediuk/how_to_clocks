
# 03 - Effects of feature selection on clock accuracy

# Shows how EWAS and aligning the genome can impact clock accuracy

#  1 - Prep workspace ====

# Load libraries
library(tidyverse)
library(cowplot)

# Load clock function
source('functions/fitClock.R')

# 2 - Load data ====

# Alignment to polar bear genome
align <- readRDS('input/pb_alignment.rds')

# Normalized betas
meth_dat <- readRDS('input/norm_betas.rds')

# Relatedness data
sibs <- readRDS('input/full_sibs.rds')

# 3 - Make list of all methyl datasets ====

# Normalized betas from ewas
ewas_betas <- list.files('output/', pattern = 'meth_ewas')
meth_list <- list()
for(i in 1:length(ewas_betas)) {
  ewas_i <- readRDS(paste0('output/', ewas_betas[i]))
  # Add to list
  meth_list[[i]] <- ewas_i
  names(meth_list)[[i]] <- substr(ewas_betas[i], 1, 13)
}

# Add normalized betas aligned to genome to list
meth_list[['meth_dat_align']] <- meth_dat %>%
  select(c(sampleId:Population, align$qname))

# Also add full set of normalized betas to list
meth_list[['meth_dat_norm']] <- meth_dat

# 4 - Fit clocks for each case (using methylation dfs in list) ====

plot_list <- list()
for(i in 1:length(meth_list)) {
  
  # Check which bears have repeat samples
  multi_samps <- meth_list[[i]] %>%
    group_by(id) %>%
    summarize(n_samps = n()) %>%
    filter(n_samps > 1) %>%
    pull(id)
  
  # Remove repeats and siblings
  distct_dat <- meth_list[[i]] %>%
    filter(! id %in% c(sibs, multi_samps))
  
  # Sample test data for feature selection clock
  train_dat <- distct_dat %>%
    group_by(Population, sex, Spec, age) %>% 
    slice_sample(n = 1)
  
  # Training data
  test_dat <- distct_dat %>%
    filter(! sampleId %in% train_dat$sampleId)
  
  # Fit clock
  fs_clock <- fitClock(train_dat, test_dat)
  
  # Add plots to list
  plot_list[[names(meth_list)[[i]]]] <- fs_clock$plot
  
}

# 5 - Plot ====

plot_grid(plotlist = plot_list, ncol = 3, labels = LETTERS[1:6], label_size = 20)

# 6 - Save ====

ggsave('fig3.tiff', plot = last_plot(), 
       device = 'tiff', path = 'figures/', dpi = 300, height = 15, width = 25, units = 'cm', bg = 'white')

