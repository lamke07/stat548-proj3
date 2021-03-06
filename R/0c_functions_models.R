library(glmnet)
library(enetLTS) # RobLL
library(zeroSum) # LZS
library(caret)
library(pROC)

fit_lasso <- function(train_x, train_y, test_x, test_y, seed, alpha = 1, sim_beta_0 = FALSE, family = "binomial", standardize = TRUE){
  set.seed(123 + seed)
  # Fit the lasso
  lambda_lasso <- cv.glmnet(x = train_x, y = train_y, alpha = alpha, family = family, standardize = standardize)$lambda.min
  glm_lasso <- glmnet(x = train_x, y = train_y, family = family, lambda = lambda_lasso, alpha = alpha, standardize = standardize)
  
  # Obtain coefficients
  glm_lasso_beta <- glm_lasso$beta
  beta_full <- t(as.vector(coef(glm_lasso)))
  colnames(beta_full) = purrr::map_chr(0:ncol(train_x), ~paste0("p_",.x))
  
  # Obtain predictions
  y_pred_lasso <- predict(glm_lasso, newx = test_x, type = "response")
  
  # Obtain evaluation metrics
  eval_metrics_lasso <- compute_metrics(y_pred = y_pred_lasso, beta_pred = glm_lasso_beta,
                                        test_y = test_y, true_beta = sim_beta_0,
                                        name = "LASSO") %>%
    bind_cols(lambda = lambda_lasso, no_zeros = sum(glm_lasso_beta == 0), beta_full)
  
  return(eval_metrics_lasso)
}

fit_enetLTS  <- function(train_x, train_y, test_x, test_y, seed, alpha = 1, sim_beta_0 = FALSE, family = "binomial", ncores = 1, standardize = TRUE){
  # Fit enetLTS
  glm_enetLTS <- enetLTS(xx = train_x, yy = as.vector(train_y), family = family, alphas = alpha, ncores = ncores, seed = 234 + seed, plot = "FALSE", scal = standardize)
  
  # Obtain coefficients
  glm_enetLTS_beta <- glm_enetLTS$coefficients
  beta_full <- t(as.vector(coefficients(glm_enetLTS)))
  colnames(beta_full) = purrr::map_chr(0:ncol(train_x), ~paste0("p_",.x))
  
  # Obtain predictions
  y_pred_LTS <- predict(glm_enetLTS, newX = test_x, type = "response", vers = "reweighted")$reweighted.response
  
  # Obtain evaluation metrics
  eval_metrics_LTS <- compute_metrics(y_pred = y_pred_LTS, beta_pred = glm_enetLTS_beta,
                                      test_y = test_y, true_beta = sim_beta_0,
                                      name = "RobLL") %>%
    bind_cols(lambda = glm_enetLTS$lambdaw, no_zeros = sum(glm_enetLTS_beta == 0), beta_full)
}

fit_zeroSum <- function(train_x, train_y, test_x, test_y, seed, alpha = 1, sim_beta_0 = FALSE, family = "binomial", standardize = TRUE){
  # LZS estimates
  set.seed(245 + seed)
  glm_zeroSum <- zeroSum(train_x, as.vector(train_y), family = family, alpha = alpha, standardize = standardize)
  
  # Obtain coefficients
  glm_zeroSum_beta <- coefficients(glm_zeroSum, s = "lambda.min")
  beta_full <- t(as.vector(glm_zeroSum_beta))
  colnames(beta_full) = purrr::map_chr(0:ncol(train_x), ~paste0("p_",.x))
  
  # Obtain predictions
  y_pred_zeroSum <- compute_mu_beta_z(beta = glm_zeroSum_beta, z = test_x)
  # y_pred_zeroSum <- predict(object = glm_zeroSum, newX = as.matrix(test_x), type = "response", s = "lambda.min")
  
  # Obtain evaluation metrics
  eval_metrics_LZS <- compute_metrics(y_pred = y_pred_zeroSum, beta_pred = glm_zeroSum_beta[-1],
                                      test_y = test_y, true_beta = sim_beta_0,
                                      name = "LZS") %>%
    bind_cols(lambda = glm_zeroSum$Lambda1SE, no_zeros = sum(beta_full[-1] == 0), beta_full)
}

fit_RobZS <- function(train_x, train_y, test_x, test_y, seed, alpha = 1, sim_beta_0 = FALSE, family = "binomial", ncores = 6, standardize = TRUE){
  # RobZS estimates
  glm_RobZS <- RobZS(xx = train_x, yy = as.vector(train_y), family = family, alphas = alpha, seed = 345 + seed, plot = FALSE, ncores = ncores, scal = standardize)
  
  # Obtain coefficients
  glm_RobZS_beta <- glm_RobZS$coefficients
  beta_full <- t(as.vector(coefficients(glm_RobZS)))
  colnames(beta_full) = purrr::map_chr(0:ncol(train_x), ~paste0("p_",.x))
  
  # Obtain predictions
  y_pred_RobZS <- predict.RobZS(glm_RobZS, newX = test_x, vers = "reweighted", type0 = "response", intercept = TRUE)$reweighted.response
  
  # Obtain evaluation metrics
  eval_metrics_RobZS <- compute_metrics(y_pred = y_pred_RobZS, beta_pred = glm_RobZS_beta,
                                        test_y = test_y, true_beta = sim_beta_0,
                                        name = "RobZS") %>%
    bind_cols(lambda = glm_RobZS$lambdaw, no_zeros = sum(glm_RobZS_beta == 0), beta_full)
}

################################################################################
################################################################################
################################################################################

train_models <- function(i, sim_train, sim_test, sim_beta_0, coeff_names, seed_select = 123, ncores = 6, standardize = TRUE){
  print(i)
  # Preprocessing
  train_x <- as.matrix(sim_train[[i]][,coeff_names])
  train_y <- as.matrix(sim_train[[i]][,"y"])
  
  test_x_1 <- as.matrix(sim_test[[i]][,c("p_0",coeff_names)])
  test_x <- as.matrix(sim_test[[i]][,coeff_names])
  test_y <- as.matrix(sim_test[[i]][,"y"])
  
  seed <- i + seed_select + length(sim_beta_0) + nrow(train_x)
  
  print("Fitting LASSO...")
  eval_metrics_lasso <- fit_lasso(train_x = train_x, train_y = train_y, 
                                  test_x = test_x, test_y = test_y,
                                  seed = seed, sim_beta_0 = sim_beta_0, standardize = standardize)
  
  print("Fitting LTS (RobLL)...")
  eval_metrics_LTS <- fit_enetLTS(train_x = train_x, train_y = train_y, 
                                  test_x = test_x, test_y = test_y,
                                  seed = seed, sim_beta_0 = sim_beta_0, standardize = standardize,
                                  ncores = ncores)
  
  print("Fitting LZS...")
  eval_metrics_LZS <- fit_zeroSum(train_x = train_x, train_y = train_y,
                                  test_x = test_x_1, test_y = test_y,
                                  seed = seed, sim_beta_0 = sim_beta_0, standardize = standardize)
  
  print("Fitting RobZS...")
  eval_metrics_RobZS <- fit_RobZS(train_x = train_x, train_y = train_y,
                                  test_x = test_x, test_y = test_y,
                                  seed = seed, sim_beta_0 = sim_beta_0, standardize = standardize,
                                  ncores = ncores)
  
  
  return(bind_rows(eval_metrics_lasso, eval_metrics_LTS, eval_metrics_LZS, eval_metrics_RobZS) %>%
           mutate(id_run = i) %>%
           relocate(id_run)
         )
}

train_models_CV <- function(train_x, train_y, test_x, test_y, sim_beta_0 = FALSE, seed = 123, ncores = 6, standardize = TRUE){
  
  print("Fitting Lasso...")
  eval_metrics_lasso <- fit_lasso(train_x = train_x, train_y = train_y,
                                  test_x = test_x, test_y = test_y,
                                  seed = seed, sim_beta_0 = sim_beta_0, standardize = standardize)

  print("Fitting LTS (RobLL)...")
  eval_metrics_LTS <- fit_enetLTS(train_x = train_x, train_y = train_y,
                                  test_x = test_x, test_y = test_y,
                                  seed = seed, sim_beta_0 = sim_beta_0, standardize = standardize,
                                  ncores = ncores)

  print("Fitting LZS...")
  eval_metrics_LZS <- fit_zeroSum(train_x = train_x, train_y = train_y,
                                  test_x = cbind(1, test_x), test_y = test_y,
                                  seed = seed, sim_beta_0 = sim_beta_0, standardize = standardize)

  print("Fitting RobZS...")
  eval_metrics_RobZS <- fit_RobZS(train_x = train_x, train_y = train_y,
                                  test_x = test_x, test_y = test_y,
                                  seed = seed, sim_beta_0 = sim_beta_0, standardize = standardize,
                                  ncores = ncores)
  
  return(bind_rows(eval_metrics_lasso, eval_metrics_LTS, eval_metrics_LZS, eval_metrics_RobZS))
}


