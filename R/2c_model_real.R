################################################################################
rm(list = ls())
source("0a_functions_utils.R")

library(cvTools)
library(glmnet)
library(enetLTS)

selbal_HIV <- readr::read_csv("data/selbal_HIV.csv") %>%
  mutate(HIV_Status = ifelse(HIV_Status == "Pos", 1, 0))

# No column with more than 90% zeros
which(colSums(selbal_HIV == 0) > 135)
# 128 HIV positive
sum(selbal_HIV$HIV_Status)
# 35% 0's
sum(selbal_HIV %>% dplyr::select(-HIV_Status) == 0)/(155*60)
################################################################################

# Define x, y and z
selbal_HIV_x <- selbal_HIV %>%
  dplyr::select(-starts_with("HIV_Status"), -starts_with("MSM")) %>%
  mutate(across(everything(), ~ifelse(.x == 0, 0.5, .x)))
selbal_HIV_y <- selbal_HIV %>%
  dplyr::select(HIV_Status)
# z <- data.frame(MSM = selbal_HIV[,61])

# Get training data
selbal_train <- bind_cols(normalize_rows(selbal_HIV_x) %>%
                            dplyr::select(-row_sum),
                          selbal_HIV_y) %>%
  arrange(HIV_Status) %>%
  mutate(row_index = row_number()) %>%
  relocate(row_index, HIV_Status)

# selbal_train_x <- selbal_train %>%
#   dplyr::select(-HIV_Status)
# selbal_train_y <- selbal_HIV_y %>%
#   dplyr::select(HIV_Status)

################################################################################

n <- nrow(selbal_train)
n0 <- nrow(selbal_train %>% filter(HIV_Status == 0))
n1 <- nrow(selbal_train %>% filter(HIV_Status == 1))
p <- ncol(selbal_HIV_x)
################################################################################
################################################################################
cv_folds <- lapply(1:20, balanced_folds, n0 = n0, n1 = n1)

################################################################################
# df <- selbal_train
# fold_select <- 1
# cv_folds <- cv_folds

train_models_all <- function(df, cv_folds, fold_select){
  cv_fold <- cv_folds[[fold_select]]
  n <- nrow(df)

  res <- purrr::map_df(1:5,
                function(i){
                  cat("\n Replication:", fold_select, "| Fold: ", i)
                  
                  df_full <- df %>%
                    left_join(cv_fold, by = "row_index") %>%
                    relocate(fold)
                  
                  train_x <- df_full %>%
                    filter(fold != i) %>%
                    dplyr::select(-c(fold, row_index, HIV_Status)) %>%
                    as.matrix()
                  
                  train_y <- df_full %>%
                    filter(fold != i) %>%
                    dplyr::pull(HIV_Status)
                  
                  test_x <- df_full %>%
                    filter(fold == i) %>%
                    dplyr::select(-c(fold, row_index, HIV_Status)) %>%
                    as.matrix()
                  
                  test_y <- df_full %>%
                    filter(fold == i) %>%
                    dplyr::pull(HIV_Status)
                  
                  res_new <- train_models_cv(train_x = train_x, train_y = train_y, test_x = test_x, test_y = test_y, seed = i*fold_select*100 + 5) %>%
                    mutate(fold = fold_select, fold_test = i)
                  
                  return(res_new)
                })

  
  return(res)
}

# x <- train_models_all(df = selbal_train, fold_select = 1, cv_folds = cv_folds)
y <- purrr::map_df(1:20, ~train_models_all(df = selbal_train, fold_select = .x, cv_folds = cv_folds))

################################################################################


train_models_cv <- function(i, train_x, train_y, test_x, test_y, sim_beta_0 = FALSE, seed_select = 123){
  # Preprocessing
  test_x_1 <- cbind(1, test_x)

  seed <- seed_select + nrow(train_x)

  # Lasso estimates
  set.seed(123 + seed)
  lambda_lasso <- cv.glmnet(x = train_x, y = train_y, alpha = 1, family = "binomial")$lambda.min
  glm_lasso <- glmnet(x = train_x, y = train_y, family = "binomial", lambda = lambda_lasso, alpha = 1)
  glm_lasso_beta <- glm_lasso$beta

  y_pred_lasso <- predict(glm_lasso, newx = test_x, type = "response")

  eval_metrics_lasso <- compute_metrics(y_pred = y_pred_lasso, beta_pred = glm_lasso_beta,
                                        test_y = test_y, true_beta = sim_beta_0,
                                        name = "LASSO") %>%
    bind_cols(lambda = lambda_lasso, beta = mean(glm_lasso_beta), no_zeros = sum(glm_lasso_beta == 0))

  # # RobLL estimates
  # glm_enetLTS <- enetLTS(xx = train_x, yy = as.vector(train_y), family = "binomial", alphas = 1, ncores = 6, seed = 234 + seed, plot = "FALSE")
  # glm_enetLTS_beta <- glm_enetLTS$coefficients
  # 
  # y_pred_LTS <- predict(glm_enetLTS, newX = test_x, type = "response", vers = "reweighted")$reweighted.response
  # 
  # eval_metrics_LTS <- compute_metrics(y_pred = y_pred_LTS, beta_pred = glm_enetLTS_beta,
  #                                     test_y = test_y, true_beta = sim_beta_0,
  #                                     name = "RobLL") %>%
  #   bind_cols(lambda = glm_enetLTS$lambdaw)

  # # LZS estimates
  # set.seed(245 + seed)
  # glm_zeroSum <- zeroSum(train_x, as.vector(train_y), family = "binomial", alpha = 1)
  # glm_zeroSum_beta <- coef(glm_zeroSum, s = "lambda.min")
  #
  # y_pred_zeroSum <- compute_mu_beta_z(beta = glm_zeroSum_beta, z = test_x_1)
  # # y_pred_zeroSum <- predict(object = glm_zeroSum, newX = as.matrix(test_x), type = "response", s = "lambda.min")
  #
  # eval_metrics_LZS <- compute_metrics(y_pred = y_pred_zeroSum, beta_pred = glm_zeroSum_beta[-1],
  #                                 test_y = test_y, true_beta = sim_beta_0,
  #                                 name = "LZS") %>%
  #   bind_cols(lambda = glm_zeroSum$Lambda1SE)
  #
  # # RobZS estimates
  # glm_RobZS <- RobZS(xx = train_x, yy = as.vector(train_y), family = "binomial", alphas = 1, seed = 345 + seed, plot = FALSE)
  # glm_RobZS_beta <- glm_RobZS$coefficients
  #
  # y_pred_RobZS <- predict.RobZS(glm_RobZS, newX = test_x, vers = "reweighted", type0 = "response", intercept = TRUE)$reweighted.response
  #
  # eval_metrics_RobZS <- compute_metrics(y_pred = y_pred_RobZS, beta_pred = glm_RobZS_beta,
  #                                 test_y = test_y, true_beta = sim_beta_0,
  #                                 name = "RobZS") %>%
  #   bind_cols(lambda = glm_RobZS$lambdaw)

  # return(bind_rows(eval_metrics_lasso, eval_metrics_LTS, eval_metrics_LZS, eval_metrics_RobZS))
  return(bind_rows(eval_metrics_lasso))
}
