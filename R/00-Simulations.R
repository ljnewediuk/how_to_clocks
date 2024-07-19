
source('functions/simCpGs.R')
source('functions/fitClock.R')

# Number of observations
n_obs = 300
n_cgs = 500

# Make n random ages between min and max
age_vec <- runif(n = n_obs, min = 0, max = 30) %>%
  sort() %>% round()

# Simulate CpGs
sim_ages <- simCpGs(n_obs = n_obs, n_cgs = n_cgs, 
                    ages = age_vec, slp = .1, err_sd = 0.8) %>%
  mutate(age = age_vec, 
         chip.ID.loc = paste0('anim', row_number()),
         tissue = rep(c('Skin', 'Blood'), length.out = n_obs)) %>%
  relocate(c(chip.ID.loc, age, tissue))

# SCENARIO 1 - Classification bias ====

# Simulates the situation in which some of the CpGs vary only in one tissue but
# not the other. Using one class to predict the other should lead to less
# accurate age estimates, while selecting a random sample of both tissues should
# produce better accuracy.

# Select a random set of CpGs to be tissue-dependent
tissue_cgs <- paste0('cg', sample(1:n_cgs, size = 50, replace = F))

# If the tissue is skin, add some noise to the proportion methylation that 
# changes the relationship between methylation and age for that site
sim_tissue <- sim_ages %>%
  mutate(noise = seq(from = 0.3, to = 0.6, length.out = n_obs),
         across(tissue_cgs, function(x) ifelse(tissue == 'Skin', x + noise, x))) %>%
  select(! noise)

# The biased training set including only skin samples
train_b <- sim_tissue %>%
  filter(tissue == 'Skin')

# The biased testing set including only blood samples
test_b <- sim_tissue %>%
  filter(tissue == 'Blood')

# Training samples randomly selected from all samples
train_unb <- sim_tissue %>%
  group_by(age) %>%
  slice_sample(n = 5)

# Testing samples remaining after selecting out the training samples
test_unb <- sim_tissue %>%
  filter(! chip.ID.loc %in% train_unb$chip.ID.loc)

# Fit biased clock and test accuracy
fitClock(train_b, test_b)$plot
fitClock(train_b, test_b)$mae
fitClock(train_b, test_b)$rsq

# Fit unbiased clock and test accuracy
fitClock(train_unb, test_unb)$plot
fitClock(train_unb, test_unb)$mae
fitClock(train_unb, test_unb)$rsq

# SCENARIO 2 - Age bias ====

# Simulates unbalanced age samples, where young samples are used to predict old
# and vice-versa

# Training samples randomly selected from all samples
train_unb <- sim_ages %>%
  group_by(age) %>%
  slice_sample(n = 5)

# Testing samples remaining after selecting out the training samples
test_unb <- sim_ages %>%
  filter(! chip.ID.loc %in% train_unb$chip.ID.loc)

# Biased training samples selecting only the first n ages up to n unbiased 
# training samples
train_b <- sim_ages %>%
  head(nrow(train_unb), 150)

test_b <- sim_ages %>%
  filter(! chip.ID.loc %in% train_b$chip.ID.loc)

# Fit biased clock and test accuracy
fitClock(train_b, test_b)$plot
fitClock(train_b, test_b)$mae
fitClock(train_b, test_b)$rsq

# Fit unbiased clock and test accuracy
fitClock(train_unb, test_unb)$plot
fitClock(train_unb, test_unb)$mae
fitClock(train_unb, test_unb)$rsq

# SCENARIO 3 - Age error ====

# Add aging error to the test samples