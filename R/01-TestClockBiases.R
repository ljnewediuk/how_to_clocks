
# 01 - Effects of sampling bias on clock performance

# Shows how age bias, sex bias, tissue bias, genetic variation, and aging error 
# can affect clock performance

# 1 - Prep workspace ====

# Load libraries
library(tidyverse)
library(cowplot)

# Load clock function
source('functions/fitClock.R')

# 2 - Load data ====

# Normalized betas and sample info
meth_dat <- readRDS('input/norm_betas.rds') 

# Relatedness data
sibs <- readRDS('input/full_sibs.rds')

# 3 - Remove repeat samples and siblings ==== 

# Check which bears have repeat samples
multi_samps <- meth_dat %>%
  group_by(id) %>%
  summarize(n_samps = n()) %>%
  filter(n_samps > 1) %>%
  pull(id)

# Distinct dataset with repeats and siblings removed
distct_dat <- meth_dat %>%
  filter(! id %in% c(sibs, multi_samps))

# Check numbers by population, age, sex, etc.
distct_dat %>%
  group_by(Population) %>%
  summarise(n())

# Pick bears from 4 populations across different parts of the Arctic (2 with 
# skin/blood and 2 with muscle) to get a good mix of tissue and geography

# Set seed for sampling
set.seed(8)

# Sample populations, sexes, ages, tissues evenly for training data
clock_dat <- distct_dat %>%
  # filter(Population %in% c('WH', 'SH', 'NB', 'LS'))
  filter(Population %in% c('WH', 'DS', 'MC', 'NB'))  %>% 
  group_by(Population, sex, Spec, age) %>% 
  slice_sample(n = 2)
  
# 4 - Unbiased clock ====

# Get training data for the unbiased clock
unb_train_dat <- clock_dat %>%
  group_by(Population, Spec, age, sex) %>%
  slice_sample(n = 1)

# Hold-out set is training data
unb_test_dat <- clock_dat %>%
  filter(! sampleId %in% unb_train_dat$sampleId)

# Plot clock
unb_clock <- fitClock(unb_train_dat, unb_test_dat)$plot

# Add label
unb_grid <- plot_grid(unb_clock, ncol = 1, labels = 'A', label_size = 20)

# 5 - Population bias ====

# Set seed for sampling
set.seed(8)

# Sample bears from WH and MC (mix of ages, sexes, tissue types) and
# use them to predict ages in NB and DS

# Training and testing data for population-bias (use one muscle and one skin &  
# blood population to predict the other muscle and skin & blood population)
# Make sure sample sizes are the same
popb_train_dat <- clock_dat %>%
  filter(Population %in% c('MC', 'WH')) %>%
  ungroup() %>%
  slice_head(n = 96)

popb_test_dat <- clock_dat %>%
  filter(Population %in% c('DS', 'NB')) %>%
  ungroup() %>%
  slice_head(n = 87)

# Training and testing data for single population (use WH to predict WH)
pops_train_dat <- clock_dat %>%
  filter(Population %in% c('WH', 'NB')) %>%
  group_by(Spec, age, sex) %>%
  slice_sample(n = 1)

pops_test_dat <- clock_dat %>%
  filter(Population %in% c('WH', 'NB') 
         & ! sampleId %in% pops_train_dat$sampleId)

# Fit both clocks
pops_clock <- fitClock(pops_train_dat, pops_test_dat, '#ff6666')$plot
popb_clock <- fitClock(popb_train_dat, popb_test_dat, '#2acaea')$plot

# Make panel plot
pop_grid <- plot_grid(pops_clock, popb_clock, 
                      ncol = 1, labels = c('B', 'G'), label_size = 20)

# 6 - Sex bias ====

# Set seed for sampling
set.seed(12)

# Sample males from all populations and use them to predict female ages

# Training and testing data for sex-bias (using females to predict males),
# making sure sample sizes are the same
sexb_train_dat <- clock_dat %>%
  filter(sex == 'F') %>% 
  ungroup() %>%
  slice_head(n = 80)

sexb_test_dat <- clock_dat %>%
  filter(sex == 'M') %>%
  ungroup() %>%
  slice_head(n = 40)

# Training and testing data for single sex (using females to predict females)
sexs_train_dat <- clock_dat %>%
  filter(sex == 'F') %>%
  group_by(Population, Spec, age) %>%
  slice_sample(n = 1)

sexs_test_dat <- clock_dat %>%
  filter(sex == 'F' & ! sampleId %in% sexs_train_dat$sampleId)

# Fit both clocks
sexs_clock <- fitClock(sexs_train_dat, sexs_test_dat, '#ff6666')$plot
sexb_clock <- fitClock(sexb_train_dat, sexb_test_dat, '#2acaea')$plot

# Make panel plot
sex_grid <- plot_grid(sexs_clock, sexb_clock, 
                      ncol = 1, labels = c('C', 'H'), label_size = 20)

# 7 - Tissue bias ====

# Set seed for sampling
set.seed(4)

# Sample skin samples and use them to predict remaining muscle and blood samples

# Training and testing data for tissue-bias (using skin to predict muscle
# and blood), making sure sample sizes are the same as the single-tissue test
tisb_train_dat <- clock_dat %>%
  filter(Spec == 'Skin')

tisb_test_dat <- clock_dat %>%
  filter(Spec %in% c('Muscle', 'Blood')) %>%
  ungroup() %>%
  slice_sample(n = 60)

# Training and testing data for single sex (using females to predict females)
tiss_train_dat <- clock_dat %>%
  filter(Spec %in% c('Skin', 'Blood')) %>%
  group_by(Population, Spec, age, sex) %>%
  slice_sample(n = 1)

tiss_test_dat <- clock_dat %>%
  filter(Spec %in% c('Skin', 'Blood') & ! sampleId %in% tiss_train_dat$sampleId)

# Fit both clocks
tiss_clock <- fitClock(tiss_train_dat, tiss_test_dat, '#ff6666')$plot
tisb_clock <- fitClock(tisb_train_dat, tisb_test_dat, '#2acaea')$plot

# Make panel plot
tis_grid <- plot_grid(tiss_clock, tisb_clock, 
                      ncol = 1, labels = c('D', 'I'), label_size = 20)

# 8 - Age bias ====

# Set seed
set.seed(14)

# Sample bears under 10 and use them to predict bears over 10

# Training and testing data for age-bias (using older bears to predict younger), 
# making sure sample sizes are the same as the single-age group test
ageb_train_dat <- clock_dat %>%
  filter(age <= 10) %>%
  ungroup() %>%
  slice_sample(n = 71)

ageb_test_dat <- clock_dat %>%
  filter(age > 10) %>%
  ungroup() %>%
  slice_sample(n = 48)

# Training and testing data for single age group including only bears over 10
ages_train_dat <- clock_dat %>%
  filter(age >= 10) %>%
  group_by(Population, Spec, sex) %>%
  slice_sample(n = 10)

ages_test_dat <- clock_dat %>%
  filter(age >= 10 & ! sampleId %in% ages_train_dat$sampleId)

# Fit both clocks
ages_clock <- fitClock(ages_train_dat, ages_test_dat, '#ff6666')$plot
ageb_clock <- fitClock(ageb_train_dat, ageb_test_dat, '#2acaea')$plot

# Make panel plot
age_grid <- plot_grid(ages_clock, ageb_clock, 
                      ncol = 1, labels = c('E', 'J'), label_size = 20)

# 9 - Aging accuracy bias ====

#  Add some error to the ages in the training samples for the unbiased clock and
#  use those samples to predict age in the remaining test samples. The error is
#  currently sampled from a normal distribution with mean 0 and sd between 0.3 
#  and 3, proportional to the age of the bear (error increases with age)

# Filter out only WH bears
age_err_dat <- clock_dat %>%
  filter(Population == 'WH')

# Sample bears by age, sex, and tissue for training set
aerr0_train_dat <- age_err_dat %>%
  group_by(age, sex, Spec) %>%
  slice_sample(n = 1)

# Add some error
aerr2_train_dat <- aerr0_train_dat %>%
  mutate(age_orig = age,
         age = ifelse(age > 1, age + round(rnorm(1, 0, age/5)), age)) %>%
  relocate(age_orig, .after = age)

# Get test data (bears not in training data)
aerr_test_dat  <- age_err_dat %>%
  filter(! sampleId %in% aerr0_train_dat$sampleId)

# Fit  both clocks
aerr0_clock <- fitClock(aerr0_train_dat, aerr_test_dat, '#ff6666')$plot
aerr2_clock <- fitClock(aerr2_train_dat, aerr_test_dat, '#2acaea')$plot

# Make panel plot
aerr_grid <- plot_grid(aerr0_clock, aerr2_clock, 
                       ncol = 1, labels = c('F', 'K'), label_size = 20)

# 10 - Save final plot ====

# Join  all panels into single figure
plot_grid(unb_grid, pop_grid, sex_grid, tis_grid, age_grid, aerr_grid, 
          nrow = 1, rel_widths = c(2, 1, 1, 1, 1, 1))

# Save
ggsave('fig2.tiff', plot = last_plot(), 
       device = 'tiff', path = 'figures/', dpi = 300, height = 15, width = 48, units = 'cm', bg = 'white')


