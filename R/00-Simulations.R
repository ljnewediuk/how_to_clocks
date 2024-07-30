
# 1 - Prep workspace ====

# Load libraries
library(tidyverse)

# Load functions
source('functions/simCpGs.R')
source('functions/fitClock.R')

# Define function to add noise to tissue and assign new df
noise_fun <- function(df, cg_sample) {
  sim_noise <- df %>%
    mutate(noise = seq(from = 0.3, to = 0.6, length.out = n_obs),
           across(cg_sample, function(x) ifelse(tissue == 'Skin', x + noise, x))) %>%
    select(! noise)
  return(sim_noise)
}

# 2 - Set up simulation variables ====

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
         id = paste0('anim', row_number()),
         tissue = rep(c('Skin', 'Blood'), length.out = n_obs)) %>%
  relocate(c(id, age, tissue))

# 3 - Bias sampling scenarios ====

# SCENARIO 1 - Classification bias

# Simulates the situation in which CpGs relationships vary between two tissues
# (blood and skin). Using one class to predict the other should lead to less
# accurate age estimates, while selecting a random sample of both tissues should
# produce better accuracy.

# Select a random set of CpGs to be tissue-dependent
tissue_cgs <- paste0('cg', sample(1:n_cgs, size = 50, replace = F))

# If the tissue is skin, add some noise to the proportion methylation that 
# changes the relationship between methylation and age for that site

# Simulated noise added to methylation data when tissue == skin
sim_tissue <- noise_fun(sim_ages, tissue_cgs)

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
  filter(! id %in% train_unb$id)

# Fit biased clock and test accuracy
fitClock(train_b, test_b, id_var = 'id')$plot
fitClock(train_b, test_b, id_var = 'id')$mae
fitClock(train_b, test_b, id_var = 'id')$rsq

# Fit unbiased clock and test accuracy
fitClock(train_unb, test_unb, id_var = 'id')$plot
fitClock(train_unb, test_unb, id_var = 'id')$mae
fitClock(train_unb, test_unb, id_var = 'id')$rsq

# SCENARIO 2 - Age bias

# Simulates unbalanced age samples, where young samples are used to predict old
# and vice-versa

# Training samples; 150 randomly selected from all samples
train_unb <- sim_ages %>%
  group_by(age) %>%
  slice_sample(n = 5)

# Testing samples remaining after selecting out the training samples
test_unb <- sim_ages %>%
  filter(! id %in% train_unb$id)

# Biased training samples selecting only the first n ages up to n unbiased 
# training samples
train_b <- sim_ages %>%
  head(nrow(train_unb), 150)

test_b <- sim_ages %>%
  filter(! id %in% train_b$id)

# Fit biased clock and test accuracy
fitClock(train_b, test_b, id_var = 'id')$plot
fitClock(train_b, test_b, id_var = 'id')$mae
fitClock(train_b, test_b, id_var = 'id')$rsq

# Fit unbiased clock and test accuracy
fitClock(train_unb, test_unb, id_var = 'id')$plot
fitClock(train_unb, test_unb, id_var = 'id')$mae
fitClock(train_unb, test_unb, id_var = 'id')$rsq

# SCENARIO 3 - Age error

# Add aging error to the training samples

# Training samples randomly selected from all samples
train_error <- sim_ages %>%
  group_by(age) %>%
  slice_sample(n = 5)

# Test samples randomly selected from all samples
test_error <- sim_ages %>%
  filter(! id %in% train_error$id)

# Iterate over different sd of error 
# Start empty data frame
error_df <- data.frame()
for(i in seq(from = 0, to = 5, by = 0.2)) {
  # Add age error
  test_error_i <- test_error %>%
    mutate(noise = rnorm(nrow(test_error), mean = 0, sd = i),
           age = round(age + noise)) %>%
    select(! noise)
  # Fit clock
  error_clock <- fitClock(train_error, test_error_i, id_var = 'id')
  # Row for df
  error_row <- data.frame(sd = i, mae = error_clock$mae, rsq = error_clock$rsq)
  # Bind rows
  error_df <- bind_rows(error_df, error_row)
}

# Plot
error_df %>%
  pivot_longer(cols = c(mae, rsq), names_to = 'metric', values_to = 'accuracy') %>%
  ggplot(aes(x = sd, y = accuracy)) +
  geom_point() + geom_smooth(method = 'loess') +
  facet_wrap(~ metric, scales = 'free') +
  theme(panel.background = element_rect(fill = 'white', colour = 'black'),
        panel.grid = element_blank(),
        plot.margin = unit(c(0.25, 0.25, 1, 1), 'cm'),
        axis.title.y = element_text(size = 15, colour = 'black', vjust = 3),
        axis.title.x = element_text(size = 15, colour = 'black', vjust = -3),
        axis.text = element_text(size = 15, colour = 'black'),
        strip.background = element_rect(fill = 'white', colour = NA),
        strip.text = element_text(size = 15, colour = 'black'))

# 4 - Feature selection scenarios ====

# Run feature selection on variables -- age and category -- then fit clock to
# see improvement. Vary the number of CpG site related to category and the
# strength of different CpG relationships with age, which will determine the 
# number of CpG sites removed and the point at which feature selection is  no
# longer beneficial

# For size 25:250, by 25, sample random CpG sites. Then, simulate different
# relationship by tissue (noise added to the sampled sites if skin). In a
# separate df, replace all the sampled CpGs with random noise not correlated
# with age. Then fit a clock for both scenarios, plus two clocks (one for each)
# starting with all the potential CpGs confounded with tissue or not correlated
# with age. This should be plotted maybe as a four-panel figure with connected
# points to show how accuracy (MAE and r) first improves as you remove the
# biased CpGs, then starts to decline when the initial set becomes small

fs_accuracy <- data.frame()
for(i in seq(from = 25, to = 475, by = 25))  {
  cg_sample <- paste0('cg', sample(1:n_cgs, size = i, replace = F))
  sim_bias <- noise_fun(sim_ages, cg_sample)
  # The biased training set including only skin samples
  train_b <- sim_bias %>%
    filter(tissue == 'Skin')
  train_b_fs <- train_b %>%
    select(! cg_sample)
  # The biased testing set including only blood samples
  test_b <- sim_bias %>%
    filter(tissue == 'Blood')
  test_b_fs <- test_b %>%
    select(! cg_sample)
  # Fit clocks
  b_clock <- fitClock(train_b, test_b, id_var = 'id')
  fs_clock <- fitClock(train_b_fs, test_b_fs, id_var = 'id')
  fs_row <- data.frame(propFeatures = i/n_cgs,
             mae = c(fs_clock$mae, b_clock$mae),
             rsq = c(fs_clock$rsq, b_clock$rsq),
             featureSelect = c('yes', 'no'))
  fs_accuracy <- bind_rows(fs_accuracy, fs_row)
}

fs_accuracy %>% 
  pivot_longer(cols = c(mae, rsq), names_to = 'metric', values_to = 'Accuracy') %>%
  mutate(featureSelect = factor(featureSelect, levels = c('no', 'yes'),
                                labels = c('All features retained', 
                                           'Feature selection')),
         metric =  factor(metric, levels = c('mae', 'rsq'), 
                                             labels = c('Median absolute error',
                                                        'Pearson\'s correlation (r)'))) %>%
  ggplot(aes(x = propFeatures, y = Accuracy)) + 
  geom_smooth(method = 'loess', 
              colour = 'red4', fill = 'red1') +
  geom_point() +
  facet_grid(vars(metric), vars(featureSelect), scales = 'free') +
  theme(panel.background = element_rect(fill = 'white', colour = 'black'),
        panel.grid = element_blank(),
        plot.margin = unit(c(0.25, 0.25, 1, 1), 'cm'),
        axis.title.y = element_text(size = 15, colour = 'black', vjust = 3),
        axis.title.x = element_text(size = 15, colour = 'black', vjust = -3),
        axis.text = element_text(size = 15, colour = 'black'),
        strip.background = element_rect(fill = 'white', colour = NA),
        strip.text = element_text(size = 15, colour = 'black')) +
  labs(x = 'Proportion group-correlated features')

ggsave('figures/sim_feature_selection.tiff', plot = last_plot(), 
       device = 'tiff', dpi = 300, height = 15, width = 17, units = 'cm')
