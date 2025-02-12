
# Normalize betas from raw idat files. This step is required to process the data
# and perform quality control steps before fitting any clocks.

# 1 Install required packages ====

# SeSaMe and minfi packages from Bioconductor (note: requires BiocManager 
# package):
bioc_packs <- c('minfi', 'sesame')
# Install if one or either is not
if(length(setdiff(bioc_packs, rownames(installed.packages()))) > 0) {
  BiocManager::install(setdiff(bioc_packs, rownames(installed.packages())))
}

# Then install developmental versions (available from 
# https://github.com/shorvath/MammalianMethylationConsortium):
dev_packs <- c('HorvathMammalMethylChip40manifest', 
               'HorvathMammalMethylChip40anno.test.unknown')
# Install if one or either is not
if(length(setdiff(dev_packs, rownames(installed.packages()))) > 0) {
  install.packages(paste0('packages/', dev_packs, '_0.2.2.tar.gz'), 
                   type = 'source', repos = NULL)
}

library(tidyverse)
library(sesame)
library(minfi)

# 2 Source function for normalizing betas ====
source('functions/normalizeBetas.R')

# Run function for each batch
for(batch_no in 1:8) {
  
  normBetas(batch_no)
  
}

