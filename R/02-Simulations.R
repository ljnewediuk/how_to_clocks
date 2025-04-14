
# Run simulations and make plots for boxes: class/age biases and feature 
# selection

# 1 - Prep workspace ====

# Load libraries
library(tidyverse)
library(cowplot)

# Load functions
source('functions/simCpGs.R')
source('functions/fitClock.R')
source('functions/makePpnPlot.R')

# Define function to add noise to one class and assign new df. This adds enough
# variation to the CpGs in the class to change the relationship between
# age and DNA methylation. The simulated example  is tissue, but it could
# represent any two-class variable.
noise_fun <- function(df, cg_sample) {
  sim_noise <- df %>%
    mutate(noise = seq(from = 0.3, to = 0.4, length.out = n_obs),
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

# Simulate CpGs with a linear relationship between age and DNA methylation and 
# age. The slope for each site's age-DNA methylation relationship is randomly 
# generated from a normal distribution with mean randomly sampled between -0.1 
# and 1 and sd 0.25. In this dataset, a separate error distribution is then 
# generated for the CpG site with mean 0 and sd 0.8 to add some realistic noise 
# to the relationship.
sim_ages <- simCpGs(n_obs = n_obs, n_cgs = n_cgs, 
                    ages = age_vec, slp = .1, err_sd = 0.5) %>%
  mutate(age = age_vec, 
         id = paste0('anim', row_number()),
         tissue = rep(c('Skin', 'Blood'), length.out = n_obs)) %>%
  relocate(c(id, age, tissue))

# Simulate CpGs with a nonlinear relationship between age and DNA methylation.
sim_ages_nonlin <- simCpGs(n_obs = n_obs, n_cgs = n_cgs, 
                    ages = age_vec, slp = 2, err_sd = 0.5, nonlin = T) %>%
  mutate(age = age_vec, 
         id = paste0('anim', row_number()),
         tissue = rep(c('Skin', 'Blood'), length.out = n_obs)) %>%
  relocate(c(id, age, tissue)) %>%
  mutate(across(everything(), ~ replace(., is.na(.), 0)))

# 3 - Bias sampling scenarios ====

# SCENARIO 1 - Classification bias

# Simulates the situation in which CpGs relationships vary between two tissues
# (blood and skin). Using one class to predict the other should lead to less
# accurate age estimates, while selecting a random sample of both tissues should
# be more accurate.

# The number of CpG sites with a class bias varies from 0 to 500. The class bias
# eliminates the relationship between DNA methylation and age for a random 
# number of CpG sites in the biased class

# If the file already exists, we can read it in
if(file.exists('output/class_bias_sim.rds')){
  class_bias_sim <- readRDS('output/class_bias_sim.rds')
} else {
  iteration <- 1
  class_bias_sim <- data.frame()
  # Otherwise we'll simulate. We repeat this process 100 times to get a bootstrap 
  # distribution of MAE and R2
  while(iteration <= 100){
    class_bias <- data.frame()
    for(i in seq(from = 0, to = 500, by = 25)) {
      # If zero biased CpG sites, the sim bias df is just the normal sim ages df
      if(i == 0) {
        sim_bias <- sim_ages
      } else {
        # Sample i CpG sites
        cg_sample <- paste0('cg', sample(1:n_cgs, size = i, replace = F))
        # Add noise to the i CpG sites for skin samples only
        sim_bias <- noise_fun(sim_ages, cg_sample)
      }
      # The biased training set including only skin samples (i-many of these skin
      # samples will have the biased age-DNA meth relationship) n = 62
      train_b <- sim_bias %>%
        filter(tissue == 'Skin') %>%
        group_by(age) %>%
        sample_n(size = 2)
      # The unbiased set including mix of biased skin and unbiased blood (n = 62)
      train_unb <- sim_bias %>%
        group_by(age, tissue) %>%
        sample_n(size = 1)
      # Testing set including only blood samples (remove any blood samples that
      # ended up in the unbiased training set
      test <- sim_bias %>%
        filter(tissue == 'Blood', ! id %in% train_unb$id)
      # Fit clocks
      b_clock <- fitClock(train_b, test, id_var = 'id', mets_only = T)
      unb_clock <- fitClock(train_unb, test, id_var = 'id', mets_only = T)
      class_row <- data.frame(iteration,
                              propFeatures = i/n_cgs,
                              mae = c(b_clock$mae, unb_clock$mae),
                              rsq = c(b_clock$rsq, unb_clock$rsq),
                              cor = c(b_clock$cor, unb_clock$cor),
                              testSamp = c('bias_only', 'mixed'))
      class_bias <- bind_rows(class_bias, class_row)
    }
    # Bind iteration data frames
    class_bias_sim <- bind_rows(class_bias_sim, class_bias)
    # Set iteration forward
    iteration <- iteration + 1
  }
  # Save for plotting
  saveRDS(class_bias_sim, 'output/class_bias_sim.rds')
}

# Summarize the class bias data by iteration
class_bias_summ <- class_bias_sim %>% 
  # Pivot mean and confidence intervals into long format
  group_by(propFeatures, testSamp) %>%
  summarize(mae_mean = mean(mae), cor_mean = mean(cor), rsq_mean = mean(rsq), 
            mae_ci = (sd(mae)/sqrt(50))*2, 
            cor_ci = (sd(cor)/sqrt(50))*2,
            rsq_ci = (sd(rsq)/sqrt(50))*2) %>%
  pivot_longer(cols = c(mae_mean, cor_mean, rsq_mean, mae_ci, cor_ci, rsq_ci), 
               names_to = c('meas', '.value'),  names_sep = '_') %>%
  # Factor the sample types and accuracy measures
  mutate(testSamp = factor(testSamp, 
                           levels = c('mixed', 'bias_only'), 
                           labels = c('Mixed', 'Biased')))

# Plot changing MAE with proportion of biased samples
mae_class_bias_plot <- class_bias_summ %>%
  filter(meas == 'mae') %>%
  prpn_plot(prpn = 'propFeatures', samp = 'testSamp',
            samp_name = 'Training sample',
            y = 'Median absolute error',
            fills = c('#a1c9ff', '#ffbe99'),
            cls <- c('#3388ff', '#ff5e00'),
            lgnd = T) +
  ylim(0,10)

# Plot changing correlation with proportion of biased samples
cor_class_bias_plot <- class_bias_summ %>%
  filter(meas == 'cor') %>%
  prpn_plot(prpn = 'propFeatures', samp = 'testSamp', 
            samp_name = 'Training sample',
            y = "Pearson's correlation", 
            fills = c('#a1c9ff', '#ffbe99'),
            cls <- c('#3388ff', '#ff5e00'),
            lgnd = F) +
  ylim(0,1)

# Plot changing R2 with proportion of biased samples
rsq_class_bias_plot <- class_bias_summ %>%
  filter(meas == 'rsq') %>%
  prpn_plot(prpn = 'propFeatures', samp = 'testSamp', 
            samp_name = 'Training sample',
            y = 'R-squared', 
            fills = c('#a1c9ff', '#ffbe99'),
            cls <- c('#3388ff', '#ff5e00'),
            lgnd = F) +
  ylim(0,1)

# Make x axis label
bias_Xlab <- ggplot() + 
  geom_text(aes(x = 0, y = 0), 
            label = 'Proportion of biased CpGs in biased class', size = 5.5) + 
  theme_void()

# Make panel plot
bias_panels <- plot_grid(rsq_class_bias_plot, mae_class_bias_plot, 
                        ncol = 2, labels = c('A', 'B'), align = 'h', 
                        rel_widths = c(0.7, 1), label_size = 20)

#  Plot with x-axis label
plot_grid(bias_panels, bias_Xlab, ncol = 1, rel_heights = c(0.9, 0.05))

# Save figure
# As tiff
ggsave('figures/sim_class_bias.tiff', plot = last_plot(), bg = 'white',
       device = 'tiff', dpi = 300, height = 10, width = 22, units = 'cm')
# As svg
ggsave('figures/sim_class_bias.svg', plot = last_plot(), bg = 'white',
       device = 'svg', dpi = 300, height = 10, width = 22, units = 'cm')

# SCENARIO 2 - Age bias

# Simulates unbalanced age samples:
#   (1) young samples aged 0-15 are used to predict old 16-30;
#   (2) old samples aged 16-30 are used to predict young 0-15
#   (3) prime-aged individuals are used to predict all ages
# and vice-versa

# Training samples; ~ 150 randomly selected from all samples and all ages
train_age <- sim_ages_nonlin %>%
  group_by(age) %>%
  slice_sample(n = 5)

# Testing samples; ~ 150 remaining after selecting out the training samples
test_age <- sim_ages_nonlin %>%
  filter(! id %in% train_age$id)

# Unbiased comparison

iteration <- 1
unb_clock_df <- data.frame()
while(iteration <= 100) {
  # Training set with all ages
unb_train <- train_age %>% 
  group_by(age) %>%
  slice_sample(n = 2)
# Testing set with all ages
unb_test <- test_age %>% 
  group_by(age) %>%
  # filter(! id %in% unb_train$id) %>%
  slice_sample(n = 4)
# Fit  the clock
unb_clock <- fitClock(unb_train, unb_test, id_var = 'id', mets_only = T)
# Bind the data frame rows
unb_clock_row <- data.frame(iteration, type = 'unbiased',
                            mae = unb_clock$mae, 
                            cor = unb_clock$cor,
                            rsq = unb_clock$rsq)
unb_clock_df <- bind_rows(unb_clock_df, unb_clock_row)
# Increase by 1 iteration
iteration <- iteration + 1
}

# Bias 1 - using young samples to predict old

iteration <- 1
b1_clock_df <- data.frame()
while(iteration <= 100) {
  # Training set with young samples
  b_train1 <- train_age %>% 
    filter(age <= 15) %>%
    group_by(age) %>%
    slice_sample(n = 4)
  # Testing set with old samples
  b_test1 <- test_age %>% 
    filter(age > 15) %>%
    group_by(age) %>%
    slice_sample(n = 10)
  # Fit  the clock
  b1_clock <- fitClock(b_train1, b_test1, id_var = 'id', mets_only = T)
  # Bind the data frame rows
  b1_clock_row <- data.frame(iteration, type = 'young-train',
                             mae = b1_clock$mae, 
                             cor = b1_clock$cor,
                             rsq = b1_clock$rsq)
  b1_clock_df <- bind_rows(b1_clock_df, b1_clock_row)
  # Increase by 1 iteration
  iteration <- iteration + 1
}

# Bias 2 - using old to predict young (just reverse testing and training from
#          previous example)

iteration <- 1
b2_clock_df <- data.frame()
while(iteration <= 100) {
  # Training set with young samples
  b_train2 <- train_age %>% 
    filter(age >= 16) %>%
    slice_sample(n = 4)
  # Testing set with old samples
  b_test2 <- test_age %>% 
    filter(age < 16) %>%
    group_by(age) %>%
    slice_sample(n = 9)
  # Fit  the clock
  b2_clock <- fitClock(b_train2, b_test2, id_var = 'id', mets_only = T)
  # Bind the data frame rows
  b2_clock_row <- data.frame(iteration, type = 'old-train',
                             mae = b2_clock$mae, 
                             cor = b2_clock$cor,
                             rsq = b2_clock$rsq)
  b2_clock_df <- bind_rows(b2_clock_df, b2_clock_row)
  # Increase by 1 iteration
  iteration <- iteration + 1
}

# Bias 3 - using prime-aged samples to predict all samples

iteration <- 1
b3_clock_df <- data.frame()
while(iteration <= 100) {
  # Prime age train
  b_train3 <- train_age %>%
    filter(age < 20 & age > 5) %>%
    slice_sample(n = 4)
  # Prime age test
  b_test3 <- test_age %>%
    group_by(age) %>%
    slice_sample(n = 8)
  # Fit the clock
  b_clock3 <- fitClock(b_test3, b_train3, id_var = 'id', mets_only = T)
  # Bind the data frame rows
  b3_clock_row <- data.frame(iteration, type = 'prime-age',
                             mae = b_clock3$mae, 
                             cor = b_clock3$cor,
                             rsq = b_clock3$rsq)
  b3_clock_df <- bind_rows(b3_clock_df, b3_clock_row)
  # Increase by 1 iteration
  iteration <- iteration + 1
}

# Summarize the clocks
age_b_clocks <- bind_rows(b1_clock_df, b2_clock_df, 
                          b3_clock_df, unb_clock_df) %>%
  group_by(type) %>%
  # Mean accuracy by group for gradient plots
  mutate(mae_accuracy = mean(mae),
         rsq_accuracy = mean(rsq))

# Plot
mae_age_b_plot <- age_b_clocks %>%
  mutate(type = factor(type, levels = c('unbiased', 'young-train',
                       'old-train', 'prime-age'))) %>%
  ggplot(aes(x = type, y = mae, fill = mae_accuracy)) +
  geom_boxplot() +
  scale_fill_gradient(low = '#3388ff', high = '#002255') +
  ylab('Median absolute error') +
  ylim(0, 8.5) +
  theme(panel.background = element_rect(fill = 'white', colour = 'black'),
        panel.grid = element_blank(),
        plot.margin = unit(c(0.25, 0.25, 0.25, 1), 'cm'),
        legend.position = 'none',
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(size = 15, colour = 'black', vjust = 3),
        axis.text.x = element_text(size = 15, colour = 'black', angle = 90),
        axis.text.y = element_text(size = 15, colour = 'black'),
        axis.title.x = element_blank(),
        strip.background = element_rect(fill = 'white', colour = NA),
        strip.text = element_text(size = 15, colour = 'black'))

rsq_age_b_plot <- age_b_clocks %>%
  mutate(type = factor(type, levels = c('unbiased', 'young-train',
                                        'old-train', 'prime-age'))) %>%
  ggplot(aes(x = type, y = rsq, fill = rsq_accuracy)) +
  geom_boxplot() +
  scale_fill_gradient(low = '#552200', high = '#ff5e00') +
  ylab('R-squared') +
  ylim(0, 1) +
  theme(panel.background = element_rect(fill = 'white', colour = 'black'),
        panel.grid = element_blank(),
        plot.margin = unit(c(0.25, 0.25, 0.25, 1), 'cm'),
        legend.position = 'none',
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(size = 15, colour = 'black', vjust = 3),
        axis.text.x = element_text(size = 15, colour = 'black', angle = 90),
        axis.text.y = element_text(size = 15, colour = 'black'),
        axis.title.x = element_blank(),
        strip.background = element_rect(fill = 'white', colour = NA),
        strip.text = element_text(size = 15, colour = 'black'))

# Make panel plot
plot_grid(rsq_age_b_plot, mae_age_b_plot, 
          ncol = 2, labels = c('A', 'B'), align = 'h', 
          label_size = 20)

# Save figure
# As tiff
ggsave('figures/sim_age_bias.tiff', plot = last_plot(), 
       device = 'tiff', dpi = 300, height = 12, width = 22, units = 'cm')
# As svg
ggsave('figures/sim_age_bias.svg', plot = last_plot(), 
       device = 'svg', dpi = 300, height = 12, width = 22, units = 'cm')

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
  print(i)
  # Repeat 100 times per level of overlap
  it <- 1
  while(it <= 100) {
    # Add age error
    test_error_i <- test_error %>%
      mutate(noise = rnorm(nrow(test_error), mean = 0, sd = i),
             age = round(age + noise)) %>%
      select(! noise)
    # Fit clock
    error_clock <- fitClock(train_error, test_error_i, id_var = 'id', 
                            mets_only = T)
    # Row for df
    error_row <- data.frame(iteration = it,
                            sd = i, 
                            mae = error_clock$mae, 
                            rsq = error_clock$rsq)
    # Bind rows
    error_df <- bind_rows(error_df, error_row)
    # Next iteration
    it <- it + 1
  }
}

# Get mean accuracy and confidence intervals
error_means <- error_df %>%
  group_by(sd) %>%
  summarize(mae_mean = mean(mae),
            rsq_mean = mean(rsq),
            mae_lower = mean(mae) - (sd(mae)/sqrt(100))*1.96,
            mae_upper = mean(mae) + (sd(mae)/sqrt(100))*1.96,
            rsq_lower = mean(rsq) - (sd(rsq)/sqrt(100))*1.96,
            rsq_upper = mean(rsq) + (sd(rsq)/sqrt(100))*1.96)

# Set custom axis scale
custom_y <- list(scale_y_continuous(limits = c(0, 5)),
                 scale_y_continuous(limits = c(0, 1)))
# Plot
mae_age_error_plot <- error_means %>%
  ggplot(aes(x = sd, y = mae_mean)) +
  geom_ribbon(fill = '#8abbff', aes(ymin = mae_lower, ymax = mae_upper)) +
  geom_line(colour = '#3388ff') +
  geom_point(colour = '#3388ff') + 
  theme(panel.background = element_rect(fill = 'white', colour = 'black'),
        panel.grid = element_blank(),
        plot.margin = unit(c(0.25, 0.25, 0, 1), 'cm'),
        axis.title.y = element_text(size = 15, colour = 'black', vjust = 3),
        axis.title.x = element_text(size = 15, colour = 'black', vjust = -3),
        axis.text = element_text(size = 15, colour = 'black'),
        strip.background = element_rect(fill = 'white', colour = NA),
        strip.text = element_text(size = 15, colour = 'black'),
        legend.position = 'none') +
  ylim(0, 5) +
  labs(x = '', y = 'Median absolute error')

rsq_age_error_plot <- error_means %>%
  ggplot(aes(x = sd, y = rsq_mean)) +
  geom_ribbon(fill = '#ffbe99', aes(ymin = rsq_lower, ymax = rsq_upper)) +
  geom_line(colour = '#ff5e00') +
  geom_point(colour = '#ff5e00') + 
  theme(panel.background = element_rect(fill = 'white', colour = 'black'),
        panel.grid = element_blank(),
        plot.margin = unit(c(0.25, 0.25, 0, 1), 'cm'),
        axis.title.y = element_text(size = 15, colour = 'black', vjust = 3),
        axis.title.x = element_text(size = 15, colour = 'black', vjust = -3),
        axis.text = element_text(size = 15, colour = 'black'),
        strip.background = element_rect(fill = 'white', colour = NA),
        strip.text = element_text(size = 15, colour = 'black'),
        legend.position = 'none') +
  ylim(0, 1) +
  labs(x = '', y = 'R-squared')

# Make x axis label
age_error_Xlab <- ggplot() + 
  geom_text(aes(x = 0, y = 0), 
            label = 'Standard deviation of age error distribution', 
            size = 5.5) + 
  theme_void()

# Make panel plot
age_error_panels <- plot_grid(rsq_age_error_plot, mae_age_error_plot, 
                         ncol = 2, labels = c('A', 'B'), align = 'h', 
                         label_size = 20)

#  Plot with x-axis label
plot_grid(age_error_panels, age_error_Xlab, ncol = 1, rel_heights = c(0.9, 0.05))

# Save figure
# tiff
ggsave('figures/sim_age_error_bias.tiff', plot = last_plot(), bg = 'white',
       device = 'tiff', dpi = 300, height = 10, width = 22, units = 'cm')
# svg
ggsave('figures/sim_age_error_bias.svg', plot = last_plot(), bg = 'white',
       device = 'svg', dpi = 300, height = 10, width = 22, units = 'cm')

# SCENARIO 4 - Age error and sample size

# Add aging error to the training samples, and compare its influence on clock
# accuracy versus sample size

# Load file if it already exists
if(file.exists('output/sim_age_error_sample_size.rds')) {
  error_df <- readRDS('output/sim_age_error_sample_size.rds')
} else {
  
  
  # Iterate over different sd of error 
  # Start empty data frame
  error_df <- data.frame()
  
  for(N in seq(from = 50, to = 1000, by = 50)) {
    
    cat('sample size = ', N)
    
    # 2 - Set up simulation variables ====
    
    # Number of observations
    n_obs = N
    n_cgs = 500
    
    # Make n random ages between min and max
    age_vec <- runif(n = n_obs, min = 0, max = 30) %>%
      sort() %>% round()
    
    # Simulate CpGs with a linear relationship between age and DNA methylation and 
    # age.
    sim_ages <- simCpGs(n_obs = n_obs, n_cgs = n_cgs, 
                        ages = age_vec, slp = .1, err_sd = 0.5) %>%
      mutate(age = age_vec, 
             id = paste0('anim', row_number()),
             tissue = rep(c('Skin', 'Blood'), length.out = n_obs)) %>%
      relocate(c(id, age, tissue))
    
    # SCENARIO 3 - Age error
    
    # Add aging error to the training samples
    
    # Training samples randomly selected from all samples
    train_error <- sim_ages %>%
      group_by(age) %>%
      slice_sample(prop = 0.5)
    
    # Test samples randomly selected from all samples
    test_error <- sim_ages %>%
      filter(! id %in% train_error$id)
    
    
    # Repeat 100 times per level of overlap/standard deviation
    for(V in seq(from = 0, to = 5, by = 0.2)) {
      it <- 1
      cat('standard deviation = ', V)
      while(it <= 100) {
        # Add age error
        test_error_i <- test_error %>%
          mutate(noise = rnorm(nrow(test_error), mean = 0, sd = V),
                 age = round(age + noise)) %>%
          select(! noise)
        # Fit clock
        error_clock <- fitClock(train_error, test_error_i, id_var = 'id', 
                                mets_only = T)
        # Row for df
        error_row <- data.frame(iteration = it,
                                N_obs = N, 
                                sd = V,
                                mae = error_clock$mae, 
                                rsq = error_clock$rsq)
        # Bind rows
        error_df <- bind_rows(error_df, error_row)
        # Next iteration
        it <- it + 1
      }
    }
  }
  
  # Save data (takes hours to run)
  saveRDS(error_df, 'output/sim_age_error_sample_size.rds')
  
}

# Get mean accuracy and confidence intervals
error_means <- error_df %>%
  group_by(N_obs, sd) %>%
  summarize(mae_mean = mean(mae),
            rsq_mean = mean(rsq),
            mae_lower = mean(mae) - (sd(mae)/sqrt(100))*1.96,
            mae_upper = mean(mae) + (sd(mae)/sqrt(100))*1.96,
            rsq_lower = mean(rsq) - (sd(rsq)/sqrt(100))*1.96,
            rsq_upper = mean(rsq) + (sd(rsq)/sqrt(100))*1.96)

# Plot
mae_err_sample_size_plot <- error_means %>%
  filter(N_obs %in% c(50, 500, 1000)) %>%
  ggplot(aes(x = sd, y = mae_mean, size = factor(N_obs))) +
  geom_ribbon(aes(ymin = mae_lower, ymax = mae_upper), 
              alpha = 0.5, colour = NA, fill = '#8abbff') +
  scale_size_manual(values = c(0.75, 1.25, 2.25), name = 'Sample size') +
  geom_line(linewidth = 0.75, colour = '#3388ff') +
  geom_point(colour = '#3388ff') + 
  theme(panel.background = element_rect(fill = 'white', colour = 'black'),
        panel.grid = element_blank(),
        plot.margin = unit(c(0.25, 0.25, 0, 1), 'cm'),
        axis.title.y = element_text(size = 15, colour = 'black', vjust = 3),
        axis.title.x = element_text(size = 15, colour = 'black', vjust = -3),
        axis.text = element_text(size = 15, colour = 'black'),
        legend.text = element_text(size = 15, colour = 'black'),
        legend.title = element_text(size = 15, colour = 'black')) +
  ylim(0, 5.5) +
  labs(x = '', y = 'Median absolute error')

rsq_err_sample_size_plot <- error_means %>%
  filter(N_obs %in% c(50, 500, 1000)) %>%
  ggplot(aes(x = sd, y = rsq_mean, size = factor(N_obs))) +
  geom_ribbon(aes(ymin = rsq_lower, ymax = rsq_upper), 
              alpha = 0.5, colour = NA, fill = '#ffbe99') +
  scale_size_manual(values = c(0.75, 1.25, 2.25), name = 'Sample size') +
  geom_line(linewidth = 0.75, colour = '#ff5e00') +
  geom_point(colour = '#ff5e00') + 
  theme(panel.background = element_rect(fill = 'white', colour = 'black'),
        panel.grid = element_blank(),
        plot.margin = unit(c(0.25, 0.25, 0, 1), 'cm'),
        axis.title.y = element_text(size = 15, colour = 'black', vjust = 3),
        axis.title.x = element_text(size = 15, colour = 'black', vjust = -3),
        axis.text = element_text(size = 15, colour = 'black'),
        legend.text = element_text(size = 15, colour = 'black'),
        legend.title = element_text(size = 15, colour = 'black')) +
  ylim(0, 1) +
  labs(x = '', y = 'R-squared')

# Make x axis label
age_error_Xlab <- ggplot() + 
  geom_text(aes(x = 0, y = 0), 
            label = 'Standard deviation of age error distribution', 
            size = 5.5) + 
  theme_void()

# Make panel plot
age_error_panels <- plot_grid(rsq_err_sample_size_plot, mae_err_sample_size_plot, 
                              ncol = 2, labels = c('A', 'B'), align = 'h', 
                              label_size = 20)

#  Plot with x-axis label
plot_grid(age_error_panels, age_error_Xlab, ncol = 1, rel_heights = c(0.9, 0.05))

# Save figure
# tiff
ggsave('figures/sim_age_error_sample_size.tiff', plot = last_plot(), bg = 'white',
       device = 'tiff', dpi = 300, height = 10, width = 22, units = 'cm')
# svg
ggsave('figures/sim_age_error_sample_size.svg', plot = last_plot(), bg = 'white',
       device = 'svg', dpi = 300, height = 10, width = 22, units = 'cm')

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

# We'll repeat this 100 times to get a bootstrapped sample of accuracy metrics
# at each proportion
fs_accuracy <- data.frame()
for(i in seq(from = 25, to = 475, by = 25))  {
  print(i)
  # Repeat 100 times per level of overlap
  it <- 1
  while(it <= 100) {
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
    b_clock <- fitClock(train_b, test_b, id_var = 'id', mets_only = T)
    fs_clock <- fitClock(train_b_fs, test_b_fs, id_var = 'id', mets_only = T)
    fs_row <- data.frame(iteration = it,
                         propFeatures = i/n_cgs,
                         mae = c(fs_clock$mae, b_clock$mae),
                         rsq = c(fs_clock$rsq, b_clock$rsq),
                         featureSelect = c('yes', 'no'))
    fs_accuracy <- bind_rows(fs_accuracy, fs_row)
    # Next iteration
    it <- it + 1
  }
}

# Get mean accuracy and confidence intervals
fs_means <- fs_accuracy %>%
  group_by(propFeatures, featureSelect) %>%
  summarize(mae_mean = mean(mae),
            rsq_mean = mean(rsq),
            mae_lower = mean(mae) - (sd(mae)/sqrt(100))*1.96,
            mae_upper = mean(mae) + (sd(mae)/sqrt(100))*1.96,
            rsq_lower = mean(rsq) - (sd(rsq)/sqrt(100))*1.96,
            rsq_upper = mean(rsq) + (sd(rsq)/sqrt(100))*1.96)

# Plot figure - median absolute error
fs_mae_plot <- fs_means %>% 
  mutate(featureSelect = factor(featureSelect, levels = c('no', 'yes'),
                                labels = c('All features retained', 
                                           'Feature selection'))) %>%
  ggplot(aes(x = propFeatures, y = mae_mean)) + 
  geom_ribbon(aes(ymin = mae_lower, ymax = mae_upper), fill = '#a1c9ff') +
  geom_line(colour = '#3388ff') +
  geom_point(colour = '#3388ff') +
  facet_wrap(~ featureSelect, scales = 'free') +
  theme(panel.background = element_rect(fill = 'white', colour = 'black'),
        panel.grid = element_blank(),
        plot.margin = unit(c(0.25, 0.25, 0, 1), 'cm'),
        axis.title.y = element_text(size = 15, colour = 'black', vjust = 3),
        axis.title.x = element_blank(),
        axis.text = element_text(size = 15, colour = 'black'),
        strip.background = element_rect(fill = 'white', colour = NA),
        strip.text = element_text(size = 15, colour = 'black'),
        legend.position = 'none') +
  ylim(0,3) + ylab('Median absolute error')

# Plot figure - R-squared
fs_rsq_plot <- fs_means %>% 
  mutate(featureSelect = factor(featureSelect, levels = c('no', 'yes'),
                                labels = c('All features retained', 
                                           'Feature selection'))) %>%
  ggplot(aes(x = propFeatures, y = rsq_mean)) + 
  geom_ribbon(aes(ymin = rsq_lower, ymax = rsq_upper), fill = '#ffbe99') +
  geom_line(colour = '#ff5e00') +
  geom_point(colour = '#ff5e00') +
  facet_wrap(~ featureSelect, scales = 'free') +
  theme(panel.background = element_rect(fill = 'white', colour = 'black'),
        panel.grid = element_blank(),
        plot.margin = unit(c(0.25, 0.25, 0, 1), 'cm'),
        axis.title.y = element_text(size = 15, colour = 'black', vjust = 3),
        axis.title.x = element_blank(),
        axis.text = element_text(size = 15, colour = 'black'),
        strip.background = element_rect(fill = 'white', colour = NA),
        strip.text = element_blank(),
        legend.position = 'none') +
  ylim(0,1) + ylab('R-squared')

# Make x axis label
fs_Xlab <- ggplot() + 
  geom_text(aes(x = 0, y = 0), 
            label = 'Proportion of biased CpGs in biased class', size = 5.5) + 
  theme_void()

# Make panel plot
fs_panels <- plot_grid(fs_mae_plot, fs_rsq_plot, 
                         ncol = 1, labels = c('A', 'B'), align = 'v', 
                         rel_heights = c(1, .8), label_size = 20)

#  Plot with x-axis label
plot_grid(fs_panels, fs_Xlab, ncol = 1, rel_heights = c(0.9, 0.05))

# Save figure
# tiff
ggsave('figures/sim_feature_selection.tiff', plot = last_plot(), bg = 'white',
       device = 'tiff', dpi = 300, height = 15, width = 17, units = 'cm')
# svg
ggsave('figures/sim_feature_selection.svg', plot = last_plot(), bg = 'white',
       device = 'svg', dpi = 300, height = 15, width = 17, units = 'cm')

