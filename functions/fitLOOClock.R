
fitLOOClock <- function(pop_dat, fs = FALSE, N = 250) {
  
  # Separate into training and testing sets (sample by pop, tissue, sex, age)
  train_dat <- pop_dat %>%
    sample_n(size = N)
    # group_by(Population, Spec, sex, age) %>%
    # slice_sample(n = 1) %>%
    # ungroup()
  
  # Perform feature selection if needed (remove sites related to sex
  # for training data only)
  if(isTRUE(fs)) {

    sex_list <- c('M', 'F')
    for(i in 1:length(sex_list)) {
      # Fit limma (can't control for population because one population is all male)
      lim_cgs <- groupLimma(dat = train_dat, variable = 'sex', value = sex_list[i],
                            control_vars <- c('Spec'), pval = 1e-3)
      # Get CpG sites correlated with age in all categories
      if(i == 1) sex_cgs <- lim_cgs
      if(i > 1) sex_cgs <- sex_cgs[sex_cgs %in% lim_cgs]
    }
    # Remove samples with feature selection
    train_dat <- train_dat %>%
      select(sampleId:Population, any_of(sex_cgs)) 
  }
  
  test_dat <- pop_dat %>%
    filter(! sampleId %in% train_dat$sampleId)
  
  # 1 - Fit LOO clock ====
  
  # Leave one individual out from training data, fit a clock with the remaining
  # population, then use it to  predict the excluded individual's age
  loo_preds <- data.frame()
  for(i in unique(train_dat$id)) {
    # Population data with excluded individual
    loo_dat <- train_dat %>%
      filter(! id == i) 
    # Methylation data and vector of ages
    loo_meth <- loo_dat %>%
      select(! sampleId:Population) %>%
      as.matrix()
    loo_ages <- loo_dat %>%
      pull(age)
    # Excluded individual
    pred_dat <- train_dat %>%
      filter(id == i)
    # Methylation data and vector of ages
    pred_meth <- pred_dat %>%
      select(! sampleId:Population) %>%
      as.matrix()
    pred_ages <- pred_dat %>%
      pull(age)
    # Fit the LOO clock
    mod.cv.loo <- cv.glmnet(x = loo_meth, y = loo_ages, alpha = 0.5, nfolds = 10)
    # Predict and prep output
    pred_row <- data.frame(age = as.numeric(pred_ages), 
                           agePredict = as.numeric(predict(mod.cv.loo, pred_meth)))
    loo_preds <- bind_rows(loo_preds, pred_row)
  }
  
  # 2 - Fit a hold-out clock ====
  
  # Fit a second hold-out clock comparing accuracy of the LOO clock to the
  # hold-out clock
  # Get training and testing methylation data
  train_meth <- train_dat %>%
    select(! sampleId:Population) %>%
    as.matrix()
  test_meth <- test_dat %>%
    select(! sampleId:Population) %>%
    as.matrix()
  # Remove extra cg sites if feature selection was performed on the fold
  if(isTRUE(fs)) {
    test_meth <- test_dat %>%
      select(any_of(sex_cgs)) %>%
      as.matrix()
  }
  # Vectors of training and testing ages
  train_ages <- train_dat %>%
    pull(age)
  test_ages <- test_dat %>%
    pull(age)
  # Fit the hold-out clock
  mod.cv.ho <- cv.glmnet(x = train_meth, y = train_ages, alpha = 0.5, nfolds = 10)
  # Predict as output
  ho_preds <- data.frame(age = as.numeric(test_ages), 
                         agePredict = as.numeric(predict(mod.cv.ho, test_meth)))
  
  # 3 - Accuracy metrics for clocks ===
  
  # Calculate mae and r2 for LOO sample
  mae_loo <- as.numeric(ie2misc::mae(loo_preds$agePredict, loo_preds$age))
  rsq_loo <- summary(lm(agePredict ~ age, data = loo_preds))$r.squared
  
  # Calculate mae and r2 for hold-out clock
  mae_ho <- as.numeric(ie2misc::mae(ho_preds$agePredict, ho_preds$age))
  rsq_ho <- summary(lm(agePredict ~ age, data = ho_preds))$r.squared
  
  # Make data frame with results
  val_df <- data.frame(Clock = c('HO', 'LOO'),
                       r2 = c(rsq_ho, rsq_loo),
                       mae = c(mae_ho, mae_loo),
                       nTrain = nrow(train_dat))
  
  return(val_df)
}
