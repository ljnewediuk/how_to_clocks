
# 04 - Stability selection clock

# Use stability selection to find sites consistently selected
# Re-do this using training and testing data, then use the stably selected CpGs
# to predict age in testing data vs. a classic LASSO model

#  1 - Prep workspace ====

# Load libraries
library(tidyverse)
library(glmnet)

# Set seed
set.seed(1234)

# Load functions
source('functions/stabSelFunctions.R')

# 2 - Load data ====

# Alignment to polar bear genome
align <- readRDS('input/pb_alignment.rds')

# Normalized betas (with sites related to population removed)
meth_dat <- readRDS('output/meth_ewas_pop.rds') %>%
  select(c(sampleId:Population, any_of(align$qname)))

# Names of cpg sites
cpg_names <- meth_dat %>%
  select(! sampleId:Population) %>% 
  colnames()

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

# Separate into training and testing sets (sample by population, tissue, etc.)
train_dat <- distct_dat %>%
  group_by(Population, Spec, sex, age) %>%
  slice_sample(n = 1)

test_dat <- distct_dat %>%
  filter(! sampleId %in% train_dat$sampleId)

# Get ages from data for feature selection (vector length = nrow(DNAm))
fs_ages <- as.numeric(train_dat$age)

# Get matrix of DNAm data
fs_meth <- train_dat %>%
  # Make chip positions rownames
  column_to_rownames('chip.ID.loc') %>%
  # Remove extra cols
  select(! sampleId:Population) %>%
  # Convert to matrix
  as.matrix()

# 4 - Stability selection ====

# Fit LASSO clock
mod.cv <- cv.glmnet(x = fs_meth, y = fs_ages, alpha = 1, nfolds = 10)

# CpGs from lasso clock (59 sites)
lasso_clock <- as.matrix(coef(mod.cv, s = 'lambda.min')) %>%
  as.data.frame() %>%
  rownames_to_column('cg') %>%
  dplyr::rename('beta' = s1) %>%
  filter(beta != 0)

# Run stability selection x1000
result <- lapply(1:1000, FUN = stabsel)

# 5 - Compute average number of selected variables (q_stab) ====

p_list <- list()
for ( i in 1 :length(result)) {
  # Get length of each result (number of selected CpGs)
  res <- as.numeric(result[[i]])
  p_list[[i]] <- length(res)
}
p_out <- do.call(cbind,p_list)
# Average number of selected CpGs
q_stab <- mean(p_out)

# 6  - Compute the selection probabilities (select_prob) ====

n_select <- rep(0, length(cpg_names))
t_select <- 0 * n_select

for (i in 1:length(result)) {
  # Get names of CpGs in result i
  tt <- names(result[[i]])[-1]
  # Add one for each selected CpG
  if(length(tt)==0){

  } else {
    t_select[which(cpg_names %in% names(result[[i]])[-1])] <- 1
    n_select <- n_select + t_select
    # Re-zero
    t_select <- 0*n_select
  }
}

# Get proportion of times each site was selected
out <- n_select/length(result)
df_out <- data.frame(cpg_names = cpg_names, pi_select = out)
# Arrange in descending order
select_prob <- df_out %>%
  arrange(desc(pi_select))

# 7 - Selection probability threshold ====

thresh <- numeric(10)
for (i in 1:10) {
  thresh[i] <- finding_thresh(q = q_stab, p = length(cpg_names), E_v = i)
}

n_cpgs <- numeric(10)
for (i in 1:10) {
  n_cpgs[i] <- nrow(select_prob[which(select_prob$pi_select > thresh[i]),])
}

# Create table of thresholds (thresh) allowing for a maximum of n false 
# discoveries (E_v) with n_cpgs
table <- data.frame(
  EV = c(1:10), 
  thresh = thresh,
  n_cpgs = n_cpgs)

# Pick clocks based on threshold error
ss_metrics <- data.frame()
for(i in 1:nrow(table)) {
  ss_clock <- ssGAM(train_dat, test_dat, select_prob, 
                    table[i ,]$thresh, by_thresh = T)
  ss_metrics_rw <- data.frame(n_cpgs = ss_clock$n_cpgs, 
                              mae = ss_clock$mae,
                              r2 = ss_clock$rsq)
  ss_metrics <- bind_rows(ss_metrics, ss_metrics_rw)
}

# Pick clocks based on n CpGs (up to 100)
ss_metrics <- data.frame()
for(i in 1:100) {
  ss_clock <- ssGAM(train_dat, test_dat, select_prob)
  ss_metrics_rw <- data.frame(n_cpgs = ss_clock$n_cpgs, 
                              mae = ss_clock$mae,
                              r2 = ss_clock$rsq)
  ss_metrics <- bind_rows(ss_metrics, ss_metrics_rw)
}

# 8 - Compare accuracy of full lasso clock vs. stability selection

# Get matrix of test DNAm data
fs_test <- test_dat %>%
  # Make chip positions rownames
  column_to_rownames('chip.ID.loc') %>%
  # Remove extra cols
  select(! sampleId:Population) %>%
  # Convert to matrix
  as.matrix()

agePredict = as.numeric(predict(mod.cv, fs_test))

# Predict using full lasso model
age_df <- data.frame(age = test_dat$age, 
                     agePredict = as.numeric(predict(mod.cv, fs_test)))
# Calculate mae 
mae <- as.numeric(ie2misc::mae(age_df$agePredict, age_df$age))
# And R2
rsq <- as.numeric(cor.test(age_df$agePredict, age_df$age)$estimate)

# Make table for lasso metrics
lasso_mets <- data.frame(n_cpgs = nrow(lasso_clock) - 1,
                         Metric = c('mae', 'r2'),
                         Accuracy = c(mae, rsq))

# Plot accuracy measures for stability selection clocks from 1-100 CpGs with
# line for accuracy of LASSO clock
ss_metrics %>%
  pivot_longer(cols = c(mae, r2), names_to = 'Metric', values_to = 'Accuracy') %>%
  ggplot(aes(x = n_cpgs, y = Accuracy)) +
  geom_smooth(method = 'loess') +
  geom_hline(data = lasso_mets, aes(yintercept = Accuracy), 
             colour = 'red', linetype = 'dashed') +
  geom_point() +
  facet_wrap(~ Metric, scales = 'free_y')

# Compare accuracy for 18 vs. 100-CpG clock (18 CpGs has better accuracy and 
# lower MAE than lasso clock with 59 CpGs)
i <- 18
ssGAM(train_dat, test_dat, select_prob)$plot
ssGAM(train_dat, test_dat, select_prob)$mae
ssGAM(train_dat, test_dat, select_prob)$rsq
i <- 100
ssGAM(train_dat, test_dat, select_prob)$plot
ssGAM(train_dat, test_dat, select_prob)$mae
ssGAM(train_dat, test_dat, select_prob)$rsq
