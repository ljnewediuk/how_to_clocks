
# After normalization, combine the methylation data with the sample info in
# preparation for making clocks.

library(tidyverse)

# 1 Load data ====

# List of matrices to load (remove first element with filtered data)
mats <- list.files('betas/', pattern = 'nbetas')

# Start character vector for batch numbers
batch_nos <- c()
# Start list for dfs
batch_list <- list()

# Set batch names (if by tissue)
# Load each matrix
for(n in 1:length(mats)) {
  mat_i <- readRDS(paste0('betas/', mats[n])) %>%
    column_to_rownames('CGid')
  # Make list with element name as matrix 1-n
  batch_list[[str_extract(mats[n], '[^.]+')]] <- mat_i
  # Make character vector of all batch numbers
  if(n %in% c(1:3, 7, 8)) b <- 1 
  if(n %in% 4:6) b <- 2
  batch_nos <- as.numeric(c(batch_nos, rep(b, length = ncol(mat_i))))
}

# 2 Arrange data and load sample information ====

# Bind betas matrices together then convert to matrix
tbetas_mat <- bind_cols(batch_list) %>%
  as.matrix()

# Make df of batch numbers and sample names from matrix rows
batch_df <- data.frame(sample = colnames(tbetas_mat), batch = batch_nos)

# Transpose betas and make df
n_betas <- tbetas_mat %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column('chip.ID.loc')

# Load updated chip sample sheets
sample_sheets <- list.files('sample_sheets/', pattern = 'updated_sample_sheet')

for(i in sample_sheets) {
  chip <- readRDS(paste0('sample_sheets/', i)) %>%
    dplyr::rename('sampleId' = Sample_Name) %>%
    dplyr::select(sampleId, chip.ID.loc)
  assign(paste0('chip', substr(i, 30, 30)), chip)
}

# Load sample information
sample_specs <- list.files('sample_sheets/', pattern = 'samples')

for(i in sample_specs) {
  assign(paste0('info', substr(i, 6, 6)), readRDS(paste0('sample_sheets/', i)))
}

# Arrange sample specs
samp_specs <- info1 %>%
  mutate(YMD = as.Date(YMD)) %>%
  bind_rows(info2, info3) %>%
  mutate(Population = 'WH') %>%
  bind_rows(info4, info5, info6, info7, info8) 

# 3 Combine array positions with sample info and methylation data ====

batch_pos <- bind_rows(chip1, chip2, chip3, chip4, 
                          chip5, chip6, chip7, chip8) %>%
  select(sampleId, chip.ID.loc) %>%
  left_join(samp_specs) %>%
  # Remove zoo bears (captive) for the new clock
  filter(! Population == 'Zoo') %>%
  # Join with methylation data
  left_join(n_betas) %>%
  # Remove duplicates
  filter(! is.na(cg00013189)) %>%
  # Fix mislabeled blood sample
  mutate(Spec = ifelse(Spec == '_Bloo', 'Blood', Spec))

# 4 Save normalized betas with info for further processing ====

saveRDS(batch_pos, 'input/norm_betas.rds')
