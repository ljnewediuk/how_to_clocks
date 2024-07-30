
library(tidyverse)
library(limma)

fitLimma <- function(meth_dat, form = c('Spec', 'age', 'sex', 'Population')) {
  
  # Combine samples into design matrix
  design_mat <- meth_dat %>%
    select(! starts_with('cg'))
  
  # Betas in rows, samples in columns
  betas <- meth_dat %>%
    select(starts_with('cg')) %>%
    as.matrix() %>%
    t()

  # Create model matrix
  mod <- model.matrix(reformulate(form, response = NULL), data = design_mat)
  
  # Fit linear model
  cb_fit <- lmFit(betas, mod)
  
  # Get components of model fit
  cb_ebayes <- eBayes(cb_fit)
  
  # Organize by pval and coefficients from lm by group
  cb_pvals_tech <- as.data.frame(cb_ebayes$p.value) %>%
    # Rename pvalue columns
    rename_with(.fn = ~paste(., 'pval', sep = '_')) %>%
    # Rename coefficient columns
    cbind(as.data.frame(cb_ebayes$coefficients)) %>%
    rename_with(.fn = ~paste(., 'coeff', sep = '_'), .cols = ! ends_with('pval')) %>%
    # Make Cg names into column
    rownames_to_column('CGid') 
  
  # Return data frame with p-values
  return(cb_pvals_tech)
  
}
