
# Clock validation procedures. Simulate the process of sampling a population and
# validating a clock using leave-one-out cross-validation versus a true hold-out 
# set. 

#  1 - Prep workspace ====

# Load libraries
library(tidyverse)
library(glmnet)
library(cowplot)

# Source functions
source('functions/fitLOOClock.R')
source('functions/fitLimma.R')
source('functions/groupLimma.R')

# 2 - Load data ====

# Alignment to polar bear genome
align <- readRDS('input/pb_alignment.rds')

# Best feature selection (no relationship with sex)
fs_sites <- readRDS('output/meth_ewas_sex.rds')

# Read in the detection p-values
dpvals <- readRDS('input/low_qc_positions.rds')

# Normalized betas (with sites removed that don't align to the genome,
# and low-qc samples)
meth_dat <- readRDS('input/norm_betas.rds') %>%
  select(sampleId:Population, any_of(align$qname)) %>%
  filter(! chip.ID.loc %in% dpvals$chip.ID.loc)

# 3 - Prep data for clock ====

# Check which bears have repeat samples
multi_samps <- meth_dat %>%
  group_by(id) %>%
  summarize(n_samps = n()) %>%
  filter(n_samps > 2) %>%
  pull(id)

# Remove repeats and siblings
distct_dat <- meth_dat %>%
  filter(! id %in% multi_samps)

# Filter feature selection sites unbiased by sex
fs_dat <- distct_dat %>%
  select(sampleId:Population, any_of(fs_sites)) 

# Take 100 samples of 400 bears from the data and make a list for comparing
# validation types, both for the sets with feature selection already applied
# and for the set without feature selection
it <- 1
# List for samples with feature selection
meth_pops_fs <- list()
# Without feature selection
meth_pops_no_fs <- list()
while(it <= 100) {
  meth_pops_fs[[it]] <- sample_n(fs_dat, size = 400)
  meth_pops_no_fs[[it]] <- distct_dat %>%
    filter(sampleId %in% meth_pops_fs[[it]]$sampleId)
  it <- it + 1
}

# Fit the hold-out and leave-one-out clocks for feature pre-selection
# Load the building file if it exists, or start a new data frame
if(file.exists('output/validation_comparison.rds')) {
  compare_vals <- readRDS('output/validation_comparison.rds')
  j <- max(compare_vals$iteration) 
} else {
  compare_vals <- data.frame()
  j <- 1
}
# Run the loop
for(i in j:length(meth_pops_fs)) {
  # Iteration
  print(i)
  # Get sample from list
  samp <- meth_pops_fs[[i]]
  # Fit LOO clock and combine data frame
  val <- fitLOOClock(samp) %>%
    mutate(iteration = i)
  compare_vals <- bind_rows(compare_vals, val)
  # Save current iteration
  saveRDS(compare_vals, 'output/validation_comparison.rds')
}

# Fit the hold-out and leave-one-out clocks without feature selection
# Load the building file if it exists, or start a new data frame
if(file.exists('output/validation_comparison_no_fs.rds')) {
  compare_vals_no_fs <- readRDS('output/validation_comparison_no_fs.rds')
  j <- max(compare_vals_no_fs$iteration) + 1
} else {
  compare_vals_no_fs <- data.frame()
  j <- 1
}
# Run the loop
for(i in j:length(meth_pops_no_fs)) {
  # Iteration
  print(i)
  # Get sample from list
  samp <- meth_pops_no_fs[[i]]
  # Fit LOO clock and combine data frame, performing feature selection on each
  # training set
  val <- fitLOOClock(samp, fs = TRUE) %>%
    mutate(iteration = i)
  compare_vals_no_fs <- bind_rows(compare_vals_no_fs, val)
  # Save current iteration
  saveRDS(compare_vals_no_fs, 'output/validation_comparison_no_fs.rds')
}

# Calculate mean accuracy by clock and accuracy metric
compare_vals_means <- compare_vals %>%
  group_by(Clock) %>%
  mutate(mean_mae = mean(mae),
         mean_r2 = mean(r2))

# 4 - Boxplots ====

# MAE
mae_plt <- compare_vals_means %>%
  ggplot(aes(x = Clock, y = mae)) +
  geom_boxplot(aes(fill = mean_mae)) +
  scale_fill_gradient(low = '#3388ff', high = '#002255') +
  ylab('Median absolute error') +
  theme(panel.background = element_rect(fill = 'white', colour = 'black'),
        panel.grid = element_blank(),
        plot.margin = unit(c(0.25, 0.25, 0.25, 1), 'cm'),
        legend.position = 'none',
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(size = 15, colour = 'black', vjust = 3),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 15, colour = 'black'),
        axis.text.x = element_text(size = 15, colour = 'black', angle = 90, hjust = 1),
        strip.background = element_rect(fill = 'white', colour = NA),
        strip.text = element_text(size = 15, colour = 'black'))

# R2
rsq_plt <- compare_vals_means %>%
  ggplot(aes(x = Clock, y = r2)) +
  geom_boxplot(aes(fill = mean_r2)) +
  scale_fill_gradient(low = '#552200', high = '#ff5e00') +
  ylab('R-squared') +
  theme(panel.background = element_rect(fill = 'white', colour = 'black'),
        panel.grid = element_blank(),
        plot.margin = unit(c(0.25, 0.25, 0.25, 1), 'cm'),
        legend.position = 'none',
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(size = 15, colour = 'black', vjust = 3),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 15, colour = 'black'),
        axis.text.x = element_text(size = 15, colour = 'black', angle = 90, hjust = 1),
        strip.background = element_rect(fill = 'white', colour = NA),
        strip.text = element_text(size = 15, colour = 'black'))

# Make panel plot
plot_grid(rsq_plt, mae_plt, 
          ncol = 2, labels = c('A', 'B'), align = 'v', label_size = 20)

# 5 - Save figure ====
# tiff
ggsave('figures/pb_validation.tiff', plot = last_plot(), 
       device = 'tiff', dpi = 300, height = 10, width = 20, units = 'cm')
# svg
ggsave('figures/pb_validation.svg', plot = last_plot(), 
       device = 'svg', dpi = 300, height = 10, width = 20, units = 'cm')
