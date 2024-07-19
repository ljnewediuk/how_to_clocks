
sampleGrps2 <- function(dat, grp_fact, grps, prop, n, complete_overlap = F) {
  
  # Get subset of data based on group factors and groups
  sub_dat <- dat %>%
    select(sampleId:Population) %>%
    filter(if_all(grp_fact, ~ . %in% grps))
  
  # Number of levels in group 2 (group 1 should have 2 levels)
  if(complete_overlap == T) {
    # Case when each group has a single level (FALSE by default)
    grp2_n <- nrow(distinct(select(sub_dat, grp_fact[2])))/2
  } else {
    grp2_n <- nrow(distinct(select(sub_dat, grp_fact[2])))
  }
  # Spread proportion over groups
  props_grps <- rep(1/grp2_n, grp2_n)
  # Vector of proportions
  props_vec <- c(props_grps*prop, rep((1-sum(props_grps*prop))/grp2_n, grp2_n))
  # Vector of number of samples
  n_samps <- round(props_vec * n)
  
  # Sample data frame
  samp_df <- sub_dat %>%
    group_by_at(grp_fact) %>% 
    nest() %>%            
    ungroup() %>% 
    arrange_at(grp_fact) %>%
    mutate(n = n_samps) %>% 
    mutate(samp = map2(data, n, sample_n)) %>% 
    select(-data) %>%
    unnest(samp) %>%
    select(! n) %>%
    left_join(dat) %>%
    relocate(grp_fact[1], .after = YMD) %>%
    relocate(grp_fact[2], .after = Born)
  
  # Return sample
  return(samp_df)
  
}
