
library(tidyverse)
library(glmnet)

fitClock <- function(train_df, test_df, id_var = 'chip.ID.loc',
                     pt_col = '#000000', pt_cat = NULL, mets_only = F) {
  # Get matrix of betas for training data
  train_m <- train_df %>%
    # Make chip positions rownames
    column_to_rownames(id_var) %>%
    # Remove extra cols
    select(starts_with('cg')) %>%
    # Convert to matrix
    as.matrix()
  
  # Get matrix of betas for test data
  test_m <- test_df %>%
    # Make chip positions rownames
    column_to_rownames(id_var) %>%
    # Remove extra cols
    select(starts_with('cg')) %>%
    # Convert to matrix
    as.matrix()
  
  # 5 Check ages ====
  
  # Add ages for training and testing
  age_df <- test_df %>%
    select(! starts_with('cg'))
  
  # Get ages for training data
  train_ages <- as.numeric(train_df$age)
  
  # 6 Fit clock and predict on training data ====
  
  # Set seed
  set.seed(321)
  
  # Glmnet model (training betas ~ ages)
  cvfit <- cv.glmnet(train_m, train_ages, nfolds = 10, alpha = .5)
  
  # Add predictions as column to ages in training data
  age_df$agePredict <- as.numeric(predict(cvfit, newx = test_m, 
                                          type = "response", s = "lambda.min"))
  
  # Unset seed
  rm(.Random.seed, envir = .GlobalEnv)
  
  # Calculate mae 
  mae <- as.numeric(ie2misc::mae(age_df$agePredict, age_df$age))
  # Correlation
  cor <- as.numeric(cor.test(age_df$agePredict, age_df$age)$estimate)
  # And R2
  rsq <- summary(lm(agePredict ~ age, data = age_df))$r.squared
  
  if(isFALSE(mets_only)) {
    # Plot
    plt <- ggplot(age_df, aes(x = age, y = agePredict)) + 
      geom_abline(intercept = 0, slope = 1, linetype = 'dashed', linewidth = 1) +
      theme(panel.background = element_rect(fill = 'white', 
                                            colour = 'black',
                                            linewidth = 1),
            plot.margin = unit(c(0.5, 0.5, 0.75, 0.75), 'cm'),
            panel.grid = element_blank(),
            axis.title.x = element_text(colour = 'black', size = 15, vjust = -3),
            axis.title.y = element_text(colour = 'black', size = 15, vjust = 5),
            axis.text = element_text(colour = 'black', size = 15),
            legend.position = 'inside',
            legend.position.inside = c(.15,.85),
            legend.title = element_blank(),
            legend.text = element_text(colour = 'black', size = 15),
            legend.key = element_rect(colour = 'black', size = 0.5)) +
      xlim(-2, 30) + ylim(-2, 30) +
      labs(x = 'Chronological age', y = 'Epigenetic age')
    
    # Add colour for plot with single colour
    if(length(pt_col) == 1) {
      col_plt <- plt +
        geom_point(size = 3, colour = pt_col) 
    }
    
    # Add colour for plot with colour by group
    if(length(pt_col) > 1) {
      col_plt <- plt +
        geom_point(size = 3, aes(colour = get(pt_cat))) +
        scale_colour_manual(values = pt_col)
    }
    
    # Combine a list with plot, age data, and accuracy measures
    age_list <- list(data = age_df, plot = col_plt, mae = mae, cor = cor, rsq = rsq)
  } else {
    age_list <- list(mae = mae, cor = cor, rsq = rsq)
  }
  return(age_list)
  
}
