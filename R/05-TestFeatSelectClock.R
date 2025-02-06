
# Effects of feature selection on clock accuracy. Cut CpG sites based on
# alignment to the polar bear genome and EWAS and then test implications for 
# polar bear clock accuracy.

#  1 - Prep workspace ====

# Load libraries
library(tidyverse)

# Load clock function
source('functions/fitClock.R')

# 2 - Load data ====

# Alignment to polar bear genome
align <- readRDS('input/pb_alignment.rds')

# Normalized betas
meth_dat <- readRDS('input/norm_betas.rds')

# 3 - Make list of all methyl datasets ====

# Normalized betas from ewas
ewas_betas <- list.files('output/', pattern = 'meth_ewas')
meth_list <- list()
for(i in 1:length(ewas_betas)) {
  ewas_i <- readRDS(paste0('output/', ewas_betas[i]))
  meth_i <- meth_dat %>%
    select(sampleId:Population, all_of(ewas_i))
  # Add to list
  meth_list[[i]] <- meth_i
  names(meth_list)[[i]] <- strsplit(ewas_betas[i], '[.]')[[1]][1]
}

# Add normalized betas aligned to genome to list
meth_list[['meth_dat_align']] <- meth_dat %>%
  select(c(sampleId:Population, align$qname))

# Also add full set of normalized betas to list
meth_list[['meth_dat_norm']] <- meth_dat

# And with both EWAS plus alignment
meth_list[['meth_dat_full_ewas']] <- meth_list[["meth_ewas_tis_sex"]] %>%
  select(c(sampleId:Population, any_of(align$qname)))

# 4 - Fit clocks for each case (using methylation dfs in list) ====

# First, get rid of bears with repeat samples and siblings

# Check which bears have repeat samples
multi_samps <- meth_dat %>%
  group_by(id) %>%
  summarize(n_samps = n()) %>%
  filter(n_samps > 1) %>%
  pull(id)
# Remove repeats and siblings
distct_dat <- meth_dat %>%
  filter(! id %in% c(sibs, multi_samps))

# Repeatedly sample bears, fit clocks, and add to list for comparison
clock_comp <- list()
it <- 1
while(it <= 500) {
  print(it)
  # Sample test data for feature selection clock
  train_samples <- distct_dat %>%
    group_by(Population, sex, Spec, age) %>% 
    slice_sample(n = 1) %>%
    pull(sampleId)
  # Training data
  test_samples <- distct_dat %>%
    filter(! sampleId %in% train_samples) %>%
    pull(sampleId)
  
  plot_list <- list()
  accuracy_list <- list()
  for(i in 1:length(meth_list)) {
    
    # Sample test data for feature selection clock
    train_dat <- meth_list[[i]] %>%
      filter(sampleId %in% train_samples)
    
    # Training data
    test_dat <- meth_list[[i]] %>%
      filter(sampleId %in% test_samples)
    
    # Fit clock
    fs_clock <- fitClock(train_df = train_dat, test_df = test_dat, mets_only = T)
    
    # Add plots to list
    plot_list[[names(meth_list)[[i]]]] <- fs_clock$plot
    accuracy_list[[names(meth_list)[[i]]]]$mae <- fs_clock$mae
    accuracy_list[[names(meth_list)[[i]]]]$rsq <- fs_clock$rsq
    
  }
  clock_comp[[it]] <- list(plot_list, accuracy_list)
  names(clock_comp[[it]]) <- c('plots', 'metrics')
  it <- it + 1
}

append_mets <- function(x) {
  x_mae_vec <- c()
  x_rsq_vec <- c()
  for(i in 1:length(x$metrics)) {
    x_mae_vec <- c(x_mae_vec, x$metrics[[i]]$mae)
    x_rsq_vec <- c(x_rsq_vec, x$metrics[[i]]$rsq)
  }
  return(
    data.frame(type = names(x$metrics),
               mae = x_mae_vec,
               rsq = x_rsq_vec)
  )
}

# Data frame with accuracies
clock_comp_df <- lapply(clock_comp, append_mets) %>%
  bind_rows()

# Get data frame of number of sites by clock type
n_sites_df <- lapply(meth_list, function(x) ncol(x)-9) %>% 
  bind_rows() %>% 
  pivot_longer(everything(), names_to = 'type', values_to = 'number_sites')

# 5 - Save ====

saveRDS(clock_comp_df, 'output/comp_fs_temp.rds')
saveRDS(n_sites_df, 'output/n_sites_fs_temp.rds')

