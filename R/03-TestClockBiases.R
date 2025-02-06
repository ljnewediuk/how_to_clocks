
# Effects of sampling bias on clock performance. Shows how age bias, sex bias, 
# tissue bias, genetic variation, and aging error can affect clock performance
# using the polar bear data.

# 1 - Prep workspace ====

# Load libraries
library(tidyverse)
library(cowplot)

# Load clock function and sampling function
source('functions/fitClock.R')
source('functions/testAccuracy.R')
source('functions/sampleGrps2.R')

# 2 - Load data ====

# Alignment to polar bear genome
align <- readRDS('input/pb_alignment.rds')

# Best feature selection (no relationship with sex)
fs_sites <- readRDS('output/meth_ewas_sex.rds')

# Read in the detection p-values
dpvals <- readRDS('input/low_qc_positions.rds')

# Normalized betas (with sites removed that don't align to the genome, 
# that are sex-dependent, and low-qc samples)
meth_dat <- readRDS('input/norm_betas.rds') %>%
  select(sampleId:Population, any_of(align$qname)) %>%
  select(sampleId:Population, any_of(fs_sites)) %>%
  filter(! chip.ID.loc %in% dpvals$chip.ID.loc)

# 3 - Population bias ====

# Filter grps of interest
pop_dat <- meth_dat %>%
  filter(Population %in% c('WH', 'SB', 'NB')) %>%
  mutate(Population = ifelse(Population %in% c('SB', 'NB'), 'BF', Population))

# Function to compare accuracy across proportion overlap between Beaufort and
# Hudson Bay population training and testing data, while drawing equal numbers 
# of male and female samples
pops_out <- testAccuracy(dat = pop_dat, 
                         grp_fact = c('Population', 'sex'),
                         grps = c('WH', 'BF', 'M', 'F'),
                         prop_list = seq(0, 1, by = 0.1),
                         train_overlap = F, test_overlap = F,
                         n_times = 100, n_train = 75, n_test = 30)

# 4 - Sex bias ====

# Function to compare accuracy across proportion overlap male and female between 
# training and testing data, while drawing equal numbers of samples from populations
sex_out <- testAccuracy(dat = meth_dat, 
                        grp_fact = c('sex', 'Population'),
                        grps = c('F', 'M', 'WH', 'SB', 'NB'),
                        prop_list = seq(0, 1, by = 0.1),
                        train_overlap = F, test_overlap = F,
                        n_times = 100, n_train = 60, n_test = 30)

# 5 - Tissue bias ====

# Filter grps of interest
tissue_dat <- meth_dat %>%
  mutate(Tissue = ifelse(Spec %in% c('Blood', 'Skin'), 'BlSk', Spec)) %>%
  relocate(Tissue, .after = Spec) %>%
  filter(! Population %in% 'NW')

# Function to compare accuracy across proportion overlap tissue between training 
# and testing data, while drawing equal numbers of samples from each population
tissue_out <- testAccuracy(dat = tissue_dat, 
                           grp_fact = c('Tissue', 'Population'),
                           grps = c('BlSk', 'Muscle', unique(tissue_dat$Population)),
                           prop_list = seq(0, 1, by = 0.1),
                           train_overlap = T, test_overlap = T,
                           n_times = 100, n_train = 100, n_test = 75)

# 6 - Age bias ====

# Filter grps of interest
age_dat <- meth_dat %>%
  mutate(ageGrp = ifelse(age <= 5, 'Immature', 'Mature')) %>%
  relocate(ageGrp, .after = age)

# Function to compare accuracy across proportion overlap immature and mature 
# bears between training and testing data, while drawing equal numbers of 
# samples from populations
age_out <- testAccuracy(dat = age_dat, 
                        grp_fact = c('ageGrp', 'Population'),
                        grps = c('Immature', 'Mature', 'WH', 'SB', 'NB'),
                        prop_list = seq(0, 1, by = 0.1),
                        train_overlap = F, test_overlap = F,
                        n_times = 100, n_train = 45, n_test = 30)


# 7 - Plot accuracy measures, add to list ====

# List of bias variables for plotting
var_list <- c('pops', 'sex', 'tissue', 'age')
# List of panel titles matching variables
panel_ts <- c(pops = 'Genetic bias', sex = 'Sex bias', 
              tissue = 'Tissue bias', age = 'Age bias')
plt_list <- list()
for(i in var_list) {
  # Get correct data
  dat <- get(paste0(i, '_out'))
  # Plot - median absolute error
  mae_plot <- dat %>% 
    group_by(overlap) %>%
    mutate(median_accuracy = median(mae)) %>%
    ggplot(aes(x = overlap, y = mae, fill = median_accuracy, group = overlap)) +
    geom_hline(yintercept = 0, linetype = 'dashed') +
    geom_boxplot() +
    scale_fill_gradient(low = '#3388ff', high = '#002255') +
    theme(panel.background = element_rect(fill = 'white', colour = 'black', linewidth = 1),
          plot.margin = unit(c(0.5, 0.5, 0, 0.75), 'cm'),
          panel.grid = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_text(colour = 'black', size = 15, vjust = 5),
          axis.text = element_text(colour = 'black', size = 15),
          strip.text = element_text(colour = 'black', size = 15),
          strip.background = element_rect(fill = 'white', colour = 'white'),
          legend.position = 'none') +
    ylab('Median absolute error')
  # Plot - R-squared
  rsq_plot <- dat %>% 
    group_by(overlap) %>%
    mutate(median_accuracy = median(r2)) %>%
    ggplot(aes(x = overlap, y = r2, fill = median_accuracy, group = overlap)) + 
    geom_hline(yintercept = 1, linetype = 'dashed') +
    geom_boxplot() +
    scale_fill_gradient(low = '#552200', high = '#ff5e00') +
    theme(panel.background = element_rect(fill = 'white', colour = 'black', linewidth = 1),
          plot.margin = unit(c(0.5, 0.5, 0, 0.75), 'cm'),
          panel.grid = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_text(colour = 'black', size = 15, vjust = 5),
          axis.text = element_text(colour = 'black', size = 15),
          strip.text = element_text(colour = 'black', size = 15),
          strip.background = element_rect(fill = 'white', colour = 'white'),
          legend.position = 'none') +
    ylim(0, 1) +
    ylab('R-squared')
    
  # Make x axis label
  bias_Xlab <- ggplot() + 
    geom_text(aes(x = 0, y = 0), 
              label = 'Proportion overlap training and testing data', size = 5.5) + 
    theme_void()
  
  # Make panel title
  bias_title <- ggplot() + 
    geom_text(aes(x = 0, y = 0), label = panel_ts[i], size = 7, fontface ='bold') + 
    theme_void()
  
  # Make panel plot
  bias_panels <- plot_grid(mae_plot, rsq_plot, 
                         ncol = 2, align = 'h', rel_widths = c(1, 1))
  
  plt <- plot_grid(bias_title,
                   bias_panels, 
                   bias_Xlab, 
                   ncol = 1, rel_heights = c(0.12, 0.9, 0.2))
  
  # Add plot to list
  plt_list[[i]] <- plt
}

# 8 - Make panel plot and save ====

# Plot panels
plot_grid(plotlist = plt_list, ncol = 2, labels = LETTERS[1:4], label_size = 20)

# Save
ggsave('pb_clock_biases.tiff', plot = last_plot(), 
       device = 'tiff', path = 'figures/', dpi = 300, height = 22, width = 45, units = 'cm', bg = 'white')


