
fitLOOClock <- function(pop_dat) {
  
  # Separate into training and testing sets (sample by pop, tissue, sex, age)
  train_dat <- pop_dat %>%
    group_by(Population, Spec, sex, age) %>%
    slice_sample(n = 1) %>%
    ungroup()
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
  
  # Sample from LOO predictions the same length as the testing data
  loo_sample <- loo_preds %>%
    sample_n(size = nrow(test_dat))
  
  # Calculate mae and r2 for LOO sample
  mae_loo <- as.numeric(ie2misc::mae(loo_sample$agePredict, loo_sample$age))
  rsq_loo <- as.numeric(cor.test(loo_sample$agePredict, loo_sample$age)$estimate)
  
  # Calculate mae and r2 for hold-out clock
  mae_ho <- as.numeric(ie2misc::mae(ho_preds$agePredict, ho_preds$age))
  rsq_ho <- as.numeric(cor.test(ho_preds$agePredict, ho_preds$age)$estimate)
  
  # Make data frame with results
  val_df <- data.frame(Clock = c('HO', 'LOO'),
                       r2 = c(rsq_ho, rsq_loo),
                       mae = c(mae_ho, mae_loo),
                       nTrain = nrow(train_dat))
  
  return(val_df)
}
