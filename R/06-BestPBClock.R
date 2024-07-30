
# Load libraries
library(tidyverse)
library(glmnet)

# Load EWAS function
source('functions/fitLimma.R')

# 2 - Load data ===

# Normalized betas and sample info
meth_dat <- readRDS('input/norm_betas.rds') %>%
  # Adjust ages to quarterly; gets closer to true ages for COYs and younger bears
  # that might be up to a year older than other individuals born in their year
  mutate(birthDate = as.Date(paste0(Born, '-01-01')),
         ageDays = as.numeric(difftime(YMD, birthDate, units = 'days')),
         correctedAge = ageDays/365,
         correctedAgeQuarterly = ifelse(correctedAge - age <= 0.25, age + 0.25, 0),
         correctedAgeQuarterly = ifelse(correctedAge - age > 0.25 & 
                                          correctedAge - age <= 0.5, age + 0.5, 
                                        correctedAgeQuarterly),
         correctedAgeQuarterly = ifelse(correctedAge - age > 0.5 & 
                                          correctedAge - age <= 0.75, age + 0.75, 
                                        correctedAgeQuarterly),
         correctedAgeQuarterly = ifelse(correctedAge - age > 0.75 & 
                                          correctedAge - age <= 0.99, age + 1, 
                                        correctedAgeQuarterly)) %>%
  select(! c(age, correctedAge, birthDate)) %>%
  rename('age' = correctedAgeQuarterly) %>%
  relocate(age, .after = id)

# 3 - Run EWAS ====

# Fit limma for each population (both sexes), each tissue (both sexes), then
# each sex (all populations and tissues) and look for CpGs correlated with age.
# Then fit a clock using just these CpGs

# Fit limma by population (controlling for age and sex if possible)
pops <- unique(meth_dat$Population)
pop_cgs <- c()
for(pop in 1:length(pops)) {
  # Get data for each population
  pop_dat <- meth_dat %>%
    filter(Population == pops[pop])
  # Fit limma, controlling also for age and sex (can't control for tissue because
  # most populations are only one type of tissue). We also can't control for sex
  # for the Norwegian Bay population because all bears are male.
  if(! pops[pop] == 'NW') {
    pop_limma <- fitLimma(meth_dat = pop_dat, form = c('age', 'sex'))
  } else {
    pop_limma <- fitLimma(meth_dat = pop_dat, form = c('age'))
  }
  # Adjusted p-values
  pvals_pop <- pop_limma %>%
    mutate(across(matches('pval'), list(adj = function(x) p.adjust(x, method = 'BH'))))
  # Get the sites correlated with age for the population
  cgs_cor_w_age <- pvals_pop %>%
    select(CGid, matches('pval_adj') & matches('age') & ! matches('Intercept')) %>%
    filter(if_any(everything(), ~ . < 0.05)) %>%
    pull(CGid)
  pop_cgs <- c(pop_cgs, cgs_cor_w_age)
  # Get only unique sites when on final population
  if(pop == length(pops)) pop_cgs <- unique(pop_cgs)
}

# Fit limma by tissue (controlling for population, age, and sex)
tissues <- unique(meth_dat$Spec)
tissue_cgs <- c()
for(tissue in 1:length(tissues)) {
  # Get data for each population
  tissue_dat <- meth_dat %>%
    filter(Spec == tissues[tissue])
  # Fit limma, controlling also for age, sex, and population
  tissue_limma <- fitLimma(meth_dat = tissue_dat, 
                           form = c('age', 'sex', 'Population'))
  # Adjusted p-values
  pvals_tissue <- tissue_limma %>%
    mutate(across(matches('pval'), list(adj = function(x) p.adjust(x, method = 'BH'))))
  # Get the sites correlated with age for the population
  cgs_cor_w_age <- pvals_tissue %>%
    select(CGid, matches('pval_adj') & matches('age') & ! matches('Intercept')) %>%
    filter(if_any(everything(), ~ . < 0.05)) %>%
    pull(CGid)
  tissue_cgs <- c(tissue_cgs, cgs_cor_w_age)
  # Get only unique sites when on final population
  if(tissue == length(tissues)) tissue_cgs <- unique(tissue_cgs)
}

# Fit limma by sex (controlling for population, age, and tissue)
sex_cgs <- c()
for(s in c('F', 'M')) {
  # Get data for each population
  sex_dat <- meth_dat %>%
    filter(sex == s)
  # Fit limma, controlling also for age, sex, and population
  sex_limma <- fitLimma(meth_dat = sex_dat, 
                           form = c('age', 'Spec', 'Population'))
  # Adjusted p-values
  pvals_sex <- sex_limma %>%
    mutate(across(matches('pval'), list(adj = function(x) p.adjust(x, method = 'BH'))))
  # Get the sites correlated with age for the population
  cgs_cor_w_age <- pvals_sex %>%
    select(CGid, matches('pval_adj') & matches('age') & ! matches('Intercept')) %>%
    filter(if_any(everything(), ~ . < 0.05)) %>%
    pull(CGid)
  sex_cgs <- c(sex_cgs, cgs_cor_w_age)
  # Get only unique sites when on final population
  if(s == 'M') sex_cgs <- unique(sex_cgs)
}

# Combine all unique age-associated CpGs
age_cgs <- c(sex_cgs, pop_cgs, tissue_cgs) %>%
  unique()

# Cg positions correlated with age
meth_dat_age <- meth_dat %>% 
  select(c(sampleId:Population, age_cgs))

# Matrix with age-correlated cg sites
meth_matrix <- meth_dat_age %>%
  select(starts_with('cg')) %>%
  as.matrix()

# Ages of samples
meth_ages <- meth_dat_age %>%
  pull(age)

# Fit clock
mod.cv <- cv.glmnet(x = meth_matrix, y = meth_ages, alpha = 0.5, nfolds = 10)

# Predict ages using model
preds <- data.frame(age = as.numeric(meth_ages), 
                       agePredict = as.numeric(predict(mod.cv, meth_matrix)))

# Get residuals for age acceleration
preds$ageAccel <- lm(preds$agePredict ~ preds$age)$resid

# Add predictions to remaining data
age_df <- meth_dat %>%
  select(! c(starts_with('cg'), starts_with('ch'), age)) %>%
  bind_cols(preds)

# Get MAE of clock
as.numeric(ie2misc::mae(age_df$agePredict, age_df$age))
# Correlation
as.numeric(cor.test(age_df$agePredict, age_df$age)$estimate)

# Plot the clock (coloured by population)
ggplot(age_df, aes(x = age, y = agePredict, colour = Population)) +
  geom_point()

# Boxplots showing variation in age acceleration among populations
age_df %>% 
  filter(Born >= 1981) %>%
  ggplot(aes(x = Population, y = ageAccel)) +
  geom_boxplot()

# Save PB ages
saveRDS(age_df, 'output/predicted_ages.rds')
