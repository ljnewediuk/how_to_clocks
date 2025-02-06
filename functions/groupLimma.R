
library(tidyverse)
library(limma)

groupLimma <- function(dat, variable, value, control_vars, pval) {
  
  # Get data for the value of interest from the variable of interest
  sub_dat <- dat %>%
    filter(get(variable) %in% value)
  # Fit limma, controlling also for control variables
  lim_mod <- fitLimma(meth_dat = sub_dat, 
                      form = c('age', control_vars))
  # Adjusted p-values
  adj_pvals <- lim_mod %>%
    mutate(across(matches('pval'), list(adj = function(x) p.adjust(x, method = 'BH'))))
  # Get the sites correlated with age for the variable
  cgs_cor_w_age <- adj_pvals %>%
    select(CGid, matches('pval_adj') & matches('age') & ! matches('Intercept')) %>%
    filter(if_any(everything(), ~ . < pval)) %>%
    pull(CGid)
  # Return CpG sites
  return(cgs_cor_w_age)
}
