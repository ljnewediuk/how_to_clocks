
testAccuracy <- function(dat, grp_fact, grps, prop_list, train_overlap, 
                         test_overlap, n_train, n_times, n_test) {
  # Get existing temporary file
  if(file.exists(paste0('temp/', grp_fact[1], '_bias_model.rds'))) {
    out <- readRDS(paste0('temp/', grp_fact[1], '_bias_model.rds'))
  } else {
    out <- data.frame()
  }
  if(nrow(out) == 0) it <- 0 else it <- max(out$iteration)
  # Repeat n times
  repeat {
    # Set number of iterations based on rep
    it <- it + 1
    print(it)
    # Load temp file if number of first iteration exceeds number of times the 
    # function should be run
    if(it > n_times) {
      warning('Number of iterations exceeds number of times the function should
            be run; returning the accuracy data')
      return(out)
    }
    # Sample and fit
    for(i in 1:length(prop_list)) {
      # Training data
      train <- sampleGrps2(dat, grp_fact, grps, 
                           prop_list[i], n = n_train, 
                           complete_overlap = train_overlap)
      # Testing data
      test <- dat %>% 
        filter(! sampleId %in% train$sampleId) %>%
        sampleGrps2(grp_fact, grps, 1, n = n_test, complete_overlap = test_overlap)
      # Fit clock
      clock <- fitClock(train, test, mets_only = T)
      # Output
      r <- data.frame(iteration = it,
                      overlap = as.numeric(prop_list[i]),
                      mae = clock$mae,
                      r2 = clock$rsq)
      out <- bind_rows(out, r)
    }
    # Save temp file with latest iteration in case the loop must be interrupted
    saveRDS(out, paste0('temp/', grp_fact[1], '_bias_model.rds'))
    # Stop loop at n_times
    if(it == n_times) break
  }
  
  return(out)
  
}

