
normBetas <- function(batch_no) {
  
  # Load sample sheet
  sample_sheet_file_name <- paste0('sample_sheets/PB_array', batch_no, '_sample_sheet', batch_no, '.rds')
  
  # Update sample sheet
  
  # List idat file names to remove any chip positions without data
  chip.IDs <- substr(list.files(paste0('iscans/batch', batch_no, '/'), recursive = T, full.names = F), 1, 12)
  chip.positions <- substr(list.files(paste0('iscans/batch', batch_no, '/'), recursive = T, full.names = F), 27, 32)
  
  # Remove any files without both red and green (i.e., two files)
  run_samples <- data.frame(chip.ID.loc = paste(chip.IDs, chip.positions, sep = '_')) %>%
    group_by(chip.ID.loc) %>%
    summarize(n_samples = n()) %>%
    filter(! n_samples == 1)
  
  # Sample sheet needs a column 'Basename' pointing to the basename of a two-colour
  # .idat file (i.e., either _Red.idat or _Grn.idat). Need to add the col Basename
  # with the file path for each sample in the corresponding row. Then use 
  # 'read.metharray.exp'to find the corresponding files using the sample sheet.
  sample_sheet <- readRDS(sample_sheet_file_name) %>%
    mutate(chip.ID.loc = paste(chip.ID, stripe, sep = '_'),
           # Add basenames (i.e., file paths for iscan files)
           Basename = paste0('iscans/batch', 
                             batch_no, '/', chip.ID, '/', chip.ID, '_', stripe)) %>%
    # Filter only iscans with both red and green iscan files
    filter(chip.ID.loc %in% run_samples$chip.ID.loc)
  
  # Create an RGChannelSet object containing raw red green channel data from .idat
  RGset <- minfi::read.metharray.exp(base = NULL, targets = sample_sheet, recursive = T) 
  
  # Annotate the RGset object with probe coordinates. 
  RGset@annotation <- c(array = 'HorvathMammalMethylChip40', annotation = "test.unknown")
  
  # Run quality control check
  
  # Detection p-value compares the methylated and unmethylated channels to
  # background signals for every site. Large detection p-values indicate 
  # poor-quality samples
  
  dpvals <- detectionP(RGset)
  
  # Get mean detection p-values
  mean_dpvals <- apply(dpvals, 2, mean) %>%
    as.data.frame() %>%
    rownames_to_column()
  
  # Make data.frame for saving and plotting
  colnames(mean_dpvals) <- c('chip.ID.loc', 'detection_p')
  
  # Make list of samples to remove
  ord_dpvals <- mean_dpvals %>%
    arrange(desc(detection_p)) 
  
  rmv_samps <- c()
  for(i in 2:nrow(ord_dpvals)-1) {
    if(ord_dpvals[i,][[2]]/ord_dpvals[i+1,][[2]] > 2) {
      rmv_samps <- c(rmv_samps, ord_dpvals[i,][[1]])
    }
  }
  
  # MethylSet object containing normalized beta values
  Mset <- minfi::preprocessNoob(RGset)
  
  # Get betas from Mset for each CG Site
  n_betas <- as_tibble(minfi::getBeta(Mset), rownames = "CGid") %>%
    # Remove the samples that didn't pass QC
    select(! all_of(rmv_samps))
  
  # Transpose normalized betas for clock-fitting
  
  # Transpose matrix after removing CGid column
  n_betas_t <- n_betas %>%
    select(! CGid) %>%
    t() %>%
    as.data.frame()
  
  # Add back CGid by renaming rest of columns to CG sites
  colnames(n_betas_t) <- n_betas$CGid
  
  # Save betas
  saveRDS(n_betas_t, paste0('betas/tbetas_PB_array', batch_no, '.rds'))
  # Save normalized betas
  saveRDS(n_betas, paste0('betas/nbetas_PB_array', batch_no, '.rds'))
  # Save updated sample sheet
  saveRDS(sample_sheet, paste0('sample_sheets/updated_sample_sheet_PB_array', batch_no, '.rds'))
}