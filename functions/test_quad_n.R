test_quad_n <- function(sdm_data, weight_resp = TRUE, 
                        weight_type = "dwpr", total_area = NULL) {
  
  # Test number of quadrature points
  train_p <- sdm_data$training$presence
  train_dat <- sdm_data$training[, !colnames(sdm_data$training) %in% "presence"]
  
  if (weight_resp) {
    if (weight_type == "normal") {
      pres <- length(train_p[train_p == 1])
      bkg <- length(train_p[train_p == 0])
      wt <- ifelse(train_p == 1, 1, pres / bkg)
      family <- "binomial"
      measure <- "auc"
    }
    if (weight_type == "iwlr") {
      wt <- (10^6)^(1 - train_p)
      family <- "binomial"
      measure <- "auc"
    }
    if (weight_type == "dwpr") {
      wt <- rep(1e-6, length(train_p))
      wt[train_p == 0] <- total_area/sum(train_p == 0)
      train_dat$presence <- train_p/wt
      family <- "poisson"
      measure <- "mse"
    }
  } else {
    wt <- NULL
    family <- "binomial"
    measure <- "auc"
  }
  
  ### Tune model
  forms <- as.formula(paste("~ 1",
                            paste(colnames(train_dat), collapse = "+"),
                            paste(paste("I(", colnames(train_dat), "^2", ")", sep = ""), 
                                  collapse = "+"), sep = "+"))
  
  training_poly <- model.matrix(forms, data = train_dat) 
  training_poly <- training_poly[,-1]
  
  lasso_cv <- glmnet::cv.glmnet(x = training_poly,
                                y = train_p,
                                family = family,
                                alpha = 1,
                                weights = wt,
                                nfolds = 5,
                                foldid = sdm_data$blocks$folds[["spatial_grid"]],
                                type.measure = measure)
  
  ### Train final model
  m_final <- glmnet::glmnet(x = training_poly,
                            y = train_p,
                            family = family,
                            alpha = 1,
                            weights = wt,
                            lambda = lasso_cv$lambda.1se)
}