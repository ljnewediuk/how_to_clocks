
simCpGs <- function(n_obs, n_cgs, ages, err_sd, nonlin = F, slp = NULL) {
  
  # Function to scale between 0 and 1
  sc01 <- function(x) {
    (x - min(x))/(max(x) - min(x))
  }
  
  # Start loop with first iteration 
  cg <- 0
  # Simulate CpG sites
  repeat {
    # Iteration
    cg <- cg + 1
    # Simulate error
    e = rnorm(n = n_obs, mean = 0, sd = err_sd)
    # Random CpGs either with relationship or not
    if(! is.null(slp)) {
      # generate a random direction, either pos or neg
      slp_dir <- sample(x = c(-slp, slp), size = 1)
      # simulate a random slope with mean == the slope
      r_slp <- sample(x = rnorm(n = 1000, mean = slp_dir, sd = 0.25), size = 1)
      # linear relationship (with error)
      y <-  sc01(ages*r_slp) + e %>%
        as.data.frame()
      # nonlinear relationship (with error)
      if(nonlin == T) {
        y <-  sc01(ages^r_slp) + e %>%
          as.data.frame()
      }
    } else {
      # no relationship
      y <- runif(n = length(ages), min = 0, max = 30) %>%
        sc01() %>% as.data.frame()
    }
    # Give the data frame column y a name
    colnames(y) <- paste0('cg', cg)
    # Bind columns
    if(cg == 1) {
      meth_array <- y
    } else {
      meth_array <- bind_cols(meth_array, y)
    }
    # End loop when reach desired number of cgs
    if(cg == n_cgs) break
  }
  # Make sure betas are capped at 1
  # meth_array <- meth_array %>%
  #   mutate(across(starts_with('cg'), function(x) ifelse(x > 1, 1, x))) %>%
  #   mutate(across(starts_with('cg'), function(x) ifelse(x < 0, 0, x)))
  
  # Return the sites
  return(meth_array)
  
}

