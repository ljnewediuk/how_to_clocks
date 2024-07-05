
# 04 - Stability selection clock

# Use stability selection to find sites consistently selected

#  1 - Prep workspace ====

# Load libraries
library(tidyverse)

# 2 - Load data ====

# Alignment to polar bear genome
align <- readRDS('input/pb_alignment.rds')

# Normalized betas (with sites related to population removed)
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

# Get ages from data for feature selection (vector length = nrow(DNAm))
fs_ages <- as.numeric(distct_dat$age)

# Get matrix of DNAm data
fs_meth <- distct_dat %>%
  # Make chip positions rownames
  column_to_rownames('chip.ID.loc') %>%
  # Remove extra cols
  select(! sampleId:Population) %>%
  # Convert to matrix
  as.matrix()

# 4 - Stability selection ====

# Fit LASSO clock
mod.cv <- cv.glmnet(x = fs_meth, y = fs_ages, alpha = 1, nfolds = 10)

# Function for stability selection
stabsel <- function(i){
  set.seed(i)
  # select subsample 
  ind_sub <- sample (1:length(fs_ages), 
                     size = floor(length(fs_ages)/2), 
                     replace = FALSE)
  fs_meth_sub <- fs_meth[ind_sub,]
  fs_ages_sub <- fs_ages[ind_sub]
  # run LASSO
  res_sub <- glmnet(x = fs_meth_sub, y = fs_ages_sub, standardize = T, alpha = 1, lambda = mod.cv$lambda.1se)
  # extract CpGs selected for prediction
  res <- coef(res_sub)
  res <- res[which(res[,1]!=0),]
  return(res)
}

result <- lapply(1:10, FUN = stabsel)
