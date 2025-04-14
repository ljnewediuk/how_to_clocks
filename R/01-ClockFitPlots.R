
# Create Fig.1 plots to demonstrate clock accuracy and the difference between a 
# clock with high R2 and high MAE vs. high R2 and low MAE.

# 1 - Prep workspace ====

# Load libraries
library(tidyverse)
library(cowplot)

# Load functions
source('functions/simCpGs.R')
source('functions/fitClock.R')

# Define function to add noise to one class and assign new df. This adds enough
# variation to the CpGs in the class to change the relationship between
# age and DNA methylation. The simulated example  is tissue, but it could
# represent any two-class variable.
noise_fun <- function(df, cg_sample, lower, upper) {
  sim_noise <- df %>%
    mutate(noise = seq(from = lower, to = upper, length.out = n_obs),
           across(cg_sample, 
                  function(x) ifelse(tissue == 'Skin', x + noise, x))) %>%
    select(! noise)
  return(sim_noise)
}

# 2 - Set up simulation variables ====

# Number of observations
n_obs = 500
n_cgs = 500

# Make n random ages between min and max
age_vec <- runif(n = n_obs, min = 0, max = 30) %>%
  sort() %>% round()

# 3 Simulate two clocks with same correlation/different MAE ====

# Both clocks have approximately the same correlation, but one clock has a
# higher MAE

# Simulate CpGs with a linear relationship between age and DNA methylation and 
# age. The slope for each site's age-DNA methylation relationship is randomly 
# generated from a normal distribution with mean randomly sampled between -0.1 
# and 1 and sd 0.25. In this dataset, a separate error distribution is then 
# generated for the CpG site with mean 0 and sd 0.8 to add some realistic noise 
# to the relationship.
sim_ages <- simCpGs(n_obs = n_obs, n_cgs = n_cgs, 
                    ages = age_vec, slp = .1, err_sd = 0.8) %>%
  mutate(age = age_vec, 
         id = paste0('anim', row_number()),
         tissue = rep(c('Skin', 'Blood'), length.out = n_obs)) %>%
  relocate(c(id, age, tissue))

cg_sample <- paste0('cg', sample(1:n_cgs, size = 475, replace = F))
# Add noise to the i CpG sites for skin samples only
sim_bias <- noise_fun(sim_ages, cg_sample, lower = 0.1, upper = 0.5)

# The biased training set which will produce a clock with high MAE
train_high_MAE <- sim_bias %>%
  filter(tissue == 'Skin') %>%
  group_by(age) %>%
  sample_n(size = 2)
# The unbiased set which will produce a clock with low MAE
train_low_MAE <- sim_bias %>%
  group_by(age, tissue) %>%
  sample_n(size = 1)
# Testing set
test <- sim_bias %>%
  filter(tissue == 'Blood', ! id %in% train_low_MAE$id)

# Fit clocks
high_MAE_clock <- fitClock(train_high_MAE, test, id_var = 'id', mets_only = F)
low_MAE_clock <- fitClock(train_low_MAE, test, id_var = 'id', mets_only = F)

# Plot panels
# High-MAE plot (will become panel C)
pC <- high_MAE_clock$plot + 
  stat_smooth(geom = 'line', lineend = 'round', 
              linewidth = 2, method = 'lm', colour = '#3388ff') +
  theme(axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())
# Low-MAE plot (will become panel D)
pD <- low_MAE_clock$plot +
  stat_smooth(geom = 'line', lineend = 'round', 
              linewidth = 2, method = 'lm', colour = '#3388ff') +
  theme(axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())

# 4 Simulate two clocks, one more accurate than the other ====

# Simulate CpGs for an accurate clock
sim_ages_accurate <- simCpGs(n_obs = n_obs, n_cgs = n_cgs, 
                     ages = age_vec, slp = .1, err_sd = 0.4) %>%
  mutate(age = age_vec, 
         id = paste0('anim', row_number()),
         tissue = rep(c('Skin', 'Blood'), length.out = n_obs)) %>%
  relocate(c(id, age, tissue))

# Accurate clock training
train_accurate <- sim_ages_accurate %>%
  group_by(age) %>%
  slice_sample(n = 10)
# Testing set for accurate clock
test_accurate <- sim_ages_accurate %>%
  filter(! id %in% train_accurate$id)

# Simulate CpGs for an inaccurate clock
sim_ages_inaccurate <- simCpGs(n_obs = n_obs, n_cgs = n_cgs, 
                             ages = age_vec, slp = .1, err_sd = 1) %>%
  mutate(age = age_vec, 
         id = paste0('anim', row_number()),
         tissue = rep(c('Skin', 'Blood'), length.out = n_obs)) %>%
  relocate(c(id, age, tissue))

# Inacurate clock training
train_inaccurate <- sim_ages_inaccurate %>%
  group_by(age) %>%
  slice_sample(n = 10)
# Testing set for accurate clock
test_inaccurate <- sim_ages_inaccurate %>%
  filter(! id %in% train_inaccurate$id)

# Fit clocks
accurate_clock <- fitClock(train_accurate, test_accurate, 
                           id_var = 'id', mets_only = F)
inaccurate_clock <- fitClock(train_inaccurate, test_inaccurate, 
                             id_var = 'id', mets_only = F)

# Plot panels
# Accurate plot (will become panel A)
pA <- accurate_clock$plot + 
  stat_smooth(geom = 'line', lineend = 'round', 
              linewidth = 2, method = 'lm', colour = '#3388ff') +
  theme(axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())
# Low-accuracy clock (will become panel B)
pB <- inaccurate_clock$plot +
  stat_smooth(geom = 'line', lineend = 'round', 
              linewidth = 2, method = 'lm', colour = '#3388ff') +
  theme(axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())

# 5 Make four-panel plot ===

# Plot four panels
p4 <- plot_grid(pA, pB, pC, pD, labels = LETTERS[1:4], label_size = 16, ncol = 2)

# X and Y labels
Ylab <- ggplot() + geom_text(aes(x = 0, y = 0), 
                             label = 'Epigenetic age', size = 7,angle = 90) + 
  theme_void()
Xlab <- ggplot() + geom_text(aes(x = 0, y = 0), 
                             label = 'Chronological age', size = 7, hjust = 0.5) + theme_void()

# Add axis labels
py <- plot_grid(Ylab, p4, rel_widths = c(0.1, 1))
plot_grid(py, Xlab, rel_heights = c(1, 0.1), ncol = 1)

# Save plot
# As tiff
ggsave('clock_intro_fig.tiff', plot = last_plot(), path = 'figures/', 
       device = 'tiff', dpi = 300, height = 12, width = 12, units = 'cm', 
       bg = 'white')

# As svg
ggsave('clock_intro_fig.svg', plot = last_plot(), path = 'figures/', 
       device = 'svg', dpi = 300, height = 18, width = 18, units = 'cm', 
       bg = 'white')
