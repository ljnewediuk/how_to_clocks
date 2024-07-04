
# 03 - Effects of feature selection on clock accuracy

# Shows how EWAS and aligning the genome can impact clock accuracy

# Load clock function
source('functions/fitClock.R')

# 2 - Load data ====

# Normalized betas and sample info
meth_dat <- readRDS('output/lm_norm_betas.rds') 

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

# Sample test data for EWAS clock
train_dat <- distct_dat %>%
  group_by(Population, sex, Spec, age) %>% 
  slice_sample(n = 1)

# Training data
test_dat <- distct_dat %>%
  filter(! sampleId %in% train_dat$sampleId)

# Fit clock
ewas_clock <- fitClock(train_dat, test_dat)

ewas_clock$plot +
  geom_jitter(aes(colour = Population))

ewas_clock$data %>%
  mutate(AgeAccel = lm(agePredict ~ age)$residuals) %>%
  ggplot(aes(x = Population, y = AgeAccel)) + geom_boxplot()

# Fit tissue-specific clocks
train_muscle <- distct_dat %>%
  filter(Spec == 'Muscle') %>%
  group_by(Population, age, sex) %>%
  slice_sample(n = 1)

test_muscle <- distct_dat %>%
  filter(Spec == 'Muscle') %>%
  filter(! sampleId %in% train_muscle$sampleId)

muscle_clock <- fitClock(train_muscle, test_muscle)

muscle_clock$plot

train_blood_skin <- distct_dat %>%
  filter(Spec %in% c('Blood', 'Skin')) %>%
  group_by(Population, age, sex) %>%
  slice_sample(n = 1)

test_blood_skin <- distct_dat %>%
  filter(Spec %in% c('Blood', 'Skin')) %>%
  filter(! sampleId %in% train_blood_skin$sampleId)

blood_skin_clock <- fitClock(train_blood_skin, test_blood_skin)

blood_skin_clock$plot

