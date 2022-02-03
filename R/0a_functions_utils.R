library(tidyverse)
library(MASS)
library(boot)
library(cvTools)

balanced_folds <- function(seed_select, n0, n1, K = 5, R = 1){
  # Returns 5 balanced folds based on grouping n0 and n1
  set.seed(seed_select + 3322)
  cv_fold0 <- cvFolds(n0, K = 5, R = 1)
  set.seed(seed_select + 2233)
  cv_fold1 <- cvFolds(n1, K = 5, R = 1)
  
  return(
    bind_rows(list(row_index = as.vector(cv_fold0$subsets), fold = cv_fold0$which),
              list(row_index = n0 + as.vector(cv_fold1$subsets), fold = cv_fold1$which)) %>%
      arrange(row_index))
}

normalize_rows <- function(df){
  # Input: data frame of containing only covariates
  stopifnot(is_tibble(df))
  
  return(
    df %>%
      mutate(row_sum = rowSums(.),
             across(!contains("row_sum"), ~ .x/row_sum))
  )
}

compute_mu_beta_z <- function(beta, z){
  # Input: vector[p] beta, matrix[n,p] z
  # Check all dimensions are the same
  stopifnot(dim(z)[2] == length(beta))
  boot::inv.logit(as.numeric(z %*% beta))
  # 1/(1+exp(-as.numeric(z %*% beta)))
}

compute_logloss <- function(y_pred, test_y){
  # Input: predictions y_pred, true values test_y
  return(-mean(test_y*log(y_pred) + (1-test_y)*log(1-y_pred)))
}

compute_metrics <- function(y_pred, beta_pred, test_y, true_beta, name, threshold = 0.5){
  y_pred_hard <- as.numeric(y_pred > threshold)
  new_error <- y_pred - test_y
  
  res_mse <- mean(new_error^2)
  res_mae <- mean(abs(new_error))
  res_ml <- compute_logloss(y_pred = y_pred, test_y = test_y)
  res_se <- caret::sensitivity(as.factor(y_pred_hard), as.factor(test_y))
  res_sp <- caret::specificity(as.factor(y_pred_hard), as.factor(test_y))
  res_auc <- as.numeric(pROC::roc(response = as.vector(test_y), predictor = as.vector(y_pred), quiet = TRUE)$auc)
  
  if(length(true_beta) == 1){
    if(true_beta == FALSE){
      return(as_tibble(list(name = name,
                            res_mse = res_mse, 
                            res_mae = res_mae,
                            res_ml = res_ml,
                            res_se = res_se,
                            res_sp = res_sp,
                            res_auc = res_auc)))
    }
  }
  
  res_fp <- sum(as.numeric(beta_pred != 0) * as.numeric(true_beta == 0))
  res_fn <- sum(as.numeric(beta_pred == 0) * as.numeric(true_beta != 0))
  
  return(as_tibble(list(name = name,
                        res_mse = res_mse, 
                        res_mae = res_mae,
                        res_ml = res_ml,
                        res_se = res_se,
                        res_sp = res_sp,
                        res_auc = res_auc,
                        res_fp = res_fp,
                        res_fn = res_fn)))
}

# compute_FP_FN <- function(beta_pred, true_beta, mode){
#   print(beta_pred)
#   print(true_beta)
#   # Compute contingency matrix for betas
#   contingency_matrix <- table(as.numeric(beta_pred != 0), as.numeric(true_beta != 0))
#   if(mode == "FP") {return(contingency_matrix[2,1])}
#   if(mode == "FN") {return(contingency_matrix[1,2])}
# }
