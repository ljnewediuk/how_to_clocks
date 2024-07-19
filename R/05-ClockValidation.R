
# 05 - Clock validation procedures

# Simulate the process of sampling a population and validating a clock using 
# leave-one-out cross-validation versus a true hold-out set. 

#  1 - Prep workspace ====

# Load libraries
library(tidyverse)
library(glmnet)

# Set seed
set.seed(1234)

# Source function
source('functions/fitLOOClock.R')

# 2 - Load data ====

# Alignment to polar bear genome
align <- readRDS('input/pb_alignment.rds')

# Normalized betas (with sites related to population removed and aligned to genome)
meth_dat <- readRDS('output/meth_ewas_pop.rds') %>%
  select(c(sampleId:Population, any_of(align$qname)))

# Relatedness data
sibs <- readRDS('input/full_sibs.rds')

# 3 - Prep data for clock ====

# Check which bears have repeat samples
multi_samps <- meth_dat %>%
  group_by(id) %>%
  summarize(n_samps = n()) %>%
  filter(n_samps > 1) %>%
  pull(id)

# Remove repeats and siblings
distct_dat <- meth_dat %>%
  filter(! id %in% c(sibs, multi_samps)) 

# Take 100 samples of 200 bears from the data and make a list for comparing
# validation types
it <- 0
meth_pops <- list()
repeat {
  it <- it + 1
  meth_pops[[it]] <- sample_n(distct_dat, size = 200)
  if(it == 100) break
}

compare_vals <- data.frame()
for(i in 1:length(meth_pops)) {
  
  print(i)
  # Get sample from list
  samp <- meth_pops[[i]]
  # Fit LOO clock
  val <- fitLOOClock(samp) %>%
    mutate(iteration = i)
  compare_vals <- bind_rows(compare_vals, val)
  
}


