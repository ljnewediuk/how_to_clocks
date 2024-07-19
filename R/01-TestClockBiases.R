
# 01 - Effects of sampling bias on clock performance

# Shows how age bias, sex bias, tissue bias, genetic variation, and aging error 
# can affect clock performance

# 1 - Prep workspace ====

# Load libraries
library(tidyverse)
library(cowplot)

# Load clock function and sampling function
source('functions/fitClock.R')
source('functions/testAccuracy.R')
source('functions/sampleGrps2.R')

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

# 4 - Population bias ====

# Filter grps of interest
pop_dat <- distct_dat %>%
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

# 5 - Sex bias ====

# Function to compare accuracy across proportion overlap male and female between 
# training and testing data, while drawing equal numbers of samples from populations
sex_out <- testAccuracy(dat = distct_dat, 
                        grp_fact = c('sex', 'Population'),
                        grps = c('F', 'M', 'WH', 'SB', 'NB'),
                        prop_list = seq(0, 1, by = 0.1),
                        train_overlap = F, test_overlap = F,
                        n_times = 100, n_train = 60, n_test = 30)

# 6 - Tissue bias ====

# Filter grps of interest
tissue_dat <- distct_dat %>%
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

# 7 - Age bias ====

# Filter grps of interest
age_dat <- distct_dat %>%
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

# 7 - Aging accuracy bias ====

#  Add some error to the ages in the training samples for the unbiased clock and
#  use those samples to predict age in the remaining test samples

# Filter training data (only distinct WH bears)
train_aerr_dat <- distct_dat %>%
  filter(Population == 'WH')

# Filter testing data (WH bears excluded from training)
test_aerr_dat <- meth_dat %>%
  filter(Population == 'WH' &
           ! sampleId %in% train_aerr_dat$sampleId) %>%
  # Remove the duplicate sample
  distinct()

# Add different degrees of error to training data and fit clock

# Define age error with sd = 5%-25% of the average maximum lifespan of the
# population (20 yrs) for mature bears >= 5 yrs
aerr_seq <- seq(.05*20, .25*20, by = 1)

# Specify number of times to repeat loop
n_times <- 100

# Fit clock for all ranges of error

# Get existing temporary file
if(file.exists('temp/age_error_model.rds')) {
  aerr_out <- readRDS('temp/age_error_model.rds')
} else {
  aerr_out <- data.frame()
}
if(nrow(aerr_out) == 0) it <- 0 else it <- max(aerr_out$iteration)
# Start loop
repeat {
  it <- it + 1
  print(it)
  # Load temp file if number of first iteration exceeds number of times the 
  # function should be run
  if(it > n_times) {
    warning('Number of iterations exceeds number of times the function should
            be run; returning the accuracy data')
    return(aerr_out)
  }
  # Fit clocks
  for(i in aerr_seq) {
    train_i <- train_aerr_dat %>%
      mutate(age = ifelse(age > 5, age + round(rnorm(1, 0, i)), age))
    clock_i <- fitClock(train_i, test_aerr_dat, '#ff6666')
    row_i <- data.frame(iteration = it, error = i, 
                        mae = clock_i$mae, r2 = clock_i$rsq)
    aerr_out <- bind_rows(aerr_out, row_i)
  }
  # Save temp file with latest iteration in case the loop must be interrupted
  saveRDS(error_table, 'temp/age_error_model.rds')
  if(it == n_times) break
}

# Fit clock without error
reg_clock <- fitClock(train_aerr_dat, test_aerr_dat, '#ff6666')

reg_clock_mets <- data.frame(metric = c('mae', 'r2'), 
                             accuracy = c(reg_clock$mae, reg_clock$rsq)) %>%
  mutate(metric_f = factor(metric, labels = c("Median absolute error", 
                                              "Pearson's correlation (r)")))

# 8 - Plot accuracy measures, add to list ====

# List of bias variables for plotting
var_list <- c('pops', 'sex', 'tissue', 'age', 'aerr')
plt_list <- list()
for(i in var_list) {
  # Get correct data
  dat <- get(paste0(i, '_out'))
  # If working with age error data, temporarily rename the "error" variable
  # for plotting
  if(i == 'aerr') dat <- rename(dat, 'overlap' = error)
  # Plot
  plt <- dat %>% 
    pivot_longer(cols = c(mae, r2), names_to = 'metric', values_to = 'accuracy') %>%
    mutate(metric_f = factor(metric, 
                             labels = c("Median absolute error", 
                                        "Pearson's correlation (r)"))) %>%
    ggplot(aes(x = factor(overlap), y = accuracy, colour = metric, fill = metric)) + 
    geom_boxplot() +
    scale_fill_manual(values = c('#0a75ad50', '#ad420a50')) +
    scale_colour_manual(values = c('#0a75ad', '#ad420a')) +
    facet_wrap(~metric_f, scales = 'free') +
    theme(panel.background = element_rect(fill = 'white', colour = 'black', linewidth = 1),
          plot.margin = unit(c(0.5, 0.5, 0.75, 0.75), 'cm'),
          panel.grid = element_blank(),
          axis.title.x = element_text(colour = 'black', size = 15, vjust = -3),
          axis.title.y = element_text(colour = 'black', size = 15, vjust = 5),
          axis.text = element_text(colour = 'black', size = 15),
          strip.text = element_text(colour = 'black', size = 15),
          strip.background = element_rect(fill = 'white', colour = 'white'),
          legend.position = 'none')
  # Add axis labels; if age error plot, add a line for the regular lasso plot
  if(i == 'aerr') {
    plt <- plt + 
      geom_hline(data = reg_clock_mets, 
                 aes(yintercept = accuracy), linetype = 'dashed') +
      labs(x = 'Standard deviation aging error (years)', y = 'Accuracy')
  } else {
    plt <- plt +
      labs(x = 'Proportion overlap training and testing data', y = 'Accuracy')
  }
  plt_list[[i]] <- plt
}

# 9 - Make panel plot and save ====

# Plot panels
plot_grid(plotlist = plt_list, ncol = 1, labels = LETTERS[1:5], label_size = 20)

# Save
ggsave('fig1.tiff', plot = last_plot(), 
       device = 'tiff', path = 'figures/', dpi = 300, height = 48, width = 24, units = 'cm', bg = 'white')


