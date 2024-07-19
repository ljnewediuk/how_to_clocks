
sampleGrps <- function(dat, grp_fact, grps, props, n) {
  
  sub_dat <- dat %>%
    select(sampleId:Population) %>%
    filter(get(grp_fact) %in% grps)
  
  n_samps <- round(props * n)
  
  samp_df <- sub_dat %>%
    group_by(get(grp_fact)) %>% 
    nest() %>%            
    ungroup() %>% 
    mutate(n = n_samps) %>% 
    mutate(samp = map2(data, n, sample_n)) %>% 
    select(-data) %>%
    unnest(samp) %>%
    select(! c(`get(grp_fact)`, n)) %>%
    left_join(dat)
  
  return(samp_df)
  
}
