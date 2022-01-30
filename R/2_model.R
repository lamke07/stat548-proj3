rm(list = ls())

setwd("~/My Drive/GitHub/RobZS/R")
source("allprogs_RobZS.R")
setwd("~/Downloads/GitHub/stat548-proj3/R")
source("0b_functions_utils.R")

library(tidyverse)
library(glmnet)
library(enetLTS) # RobLL
library(zeroSum) # LZS
library(caret)
library(pROC)

################################################################################
sim1 <- readRDS("data/sim1.RDS")

sim1_train <- do.call("rbind", lapply(sim1, '[', "Z_contaminated"))
sim1_beta <- do.call("rbind", lapply(sim1, '[', "beta"))[[1]]
sim1_test <- do.call("rbind", lapply(sim1, '[', "Z_test"))

sim1_beta_0 <- sim1_beta[-1]
coeff_names1 <- purrr::map_chr(1:length(sim1_beta_0), ~paste0("p_",.x))

################################################################################

sim2 <- readRDS("data/sim2.RDS")

sim2_train <- do.call("rbind", lapply(sim2, '[', "Z_contaminated"))
sim2_beta <- do.call("rbind", lapply(sim2, '[', "beta"))[[1]]
sim2_test <- do.call("rbind", lapply(sim2, '[', "Z_test"))

sim2_beta_0 <- sim2_beta[-1]
coeff_names2 <- purrr::map_chr(1:length(sim2_beta_0), ~paste0("p_",.x))


################################################################################
################################################################################

train_models <- function(i, sim_train, sim_test, sim_beta_0, coeff_names, seed_select = 123){
  print(i)
  # Preprocessing
  train_x <- as.matrix(sim_train[[i]][,coeff_names])
  train_y <- as.matrix(sim_train[[i]][,"y"])
  
  test_x_1 <- as.matrix(sim_test[[i]][,c("p_0",coeff_names)])
  test_x <- as.matrix(sim_test[[i]][,coeff_names])
  test_y <- as.matrix(sim_test[[i]][,"y"])
  
  seed <- i + seed_select + length(sim_beta_0) + nrow(train_x)
  
  # Lasso estimates
  set.seed(123 + seed)
  lambda_lasso <- cv.glmnet(x = train_x, y = train_y)$lambda.1se
  glm_lasso <- glmnet(x = train_x, y = train_y, family = "binomial", lambda = lambda_lasso, alpha = 1)
  glm_lasso_beta <- glm_lasso$beta
  
  y_pred_lasso <- predict(glm_lasso, newx = test_x, type = "response")
  
  eval_metrics_lasso <- compute_metrics(y_pred = y_pred_lasso, beta_pred = glm_lasso_beta,
                                  test_y = test_y, true_beta = sim_beta_0,
                                  name = "LASSO") %>%
    bind_cols(lambda = lambda_lasso)
  
  # RobLL estimates
  glm_enetLTS <- enetLTS(xx = train_x, yy = as.vector(train_y), family = "binomial", alphas = 1, ncores = 6, seed = 234 + seed, plot = "FALSE")
  glm_enetLTS_beta <- glm_enetLTS$coefficients
  
  y_pred_LTS <- predict(glm_enetLTS, newX = test_x, type = "response", vers = "reweighted")$reweighted.response
  
  eval_metrics_LTS <- compute_metrics(y_pred = y_pred_LTS, beta_pred = glm_enetLTS_beta,
                                  test_y = test_y, true_beta = sim_beta_0,
                                  name = "RobLL") %>%
    bind_cols(lambda = glm_enetLTS$lambdaw)
  
  # LZS estimates
  set.seed(245 + seed)
  glm_zeroSum <- zeroSum(train_x, as.vector(train_y), family = "binomial", alpha = 1)
  glm_zeroSum_beta <- coef(glm_zeroSum, s = "lambda.1se")
  
  y_pred_zeroSum <- compute_mu_beta_z(beta = glm_zeroSum_beta, z = test_x_1)
  # y_pred_zeroSum <- predict(object = glm_zeroSum, newX = as.matrix(test_x), type = "response", s = "lambda.1se")
  
  eval_metrics_LZS <- compute_metrics(y_pred = y_pred_zeroSum, beta_pred = glm_zeroSum_beta[-1],
                                  test_y = test_y, true_beta = sim_beta_0,
                                  name = "LZS") %>%
    bind_cols(lambda = glm_zeroSum$Lambda1SE)
  
  # RobZS estimates
  glm_RobZS <- RobZS(xx = train_x, yy = as.vector(train_y), family = "binomial", alphas = 1, seed = 345 + seed, plot = FALSE)
  glm_RobZS_beta <- glm_RobZS$coefficients
  
  y_pred_RobZS <- predict.RobZS(glm_RobZS, newX = test_x, vers = "reweighted", type0 = "response", intercept = TRUE)$reweighted.response
  
  eval_metrics_RobZS <- compute_metrics(y_pred = y_pred_RobZS, beta_pred = glm_RobZS_beta,
                                  test_y = test_y, true_beta = sim_beta_0,
                                  name = "RobZS") %>%
    bind_cols(lambda = glm_RobZS$lambdaw)
  
  return(bind_rows(eval_metrics_lasso, eval_metrics_LTS, eval_metrics_LZS, eval_metrics_RobZS))
}

################################################################################


# x1 <- purrr::map_df(1:100, ~train_models(.x, sim_train = sim1_train, sim_test = sim1_test, sim_beta_0 = sim1_beta_0, coeff_names = coeff_names1, seed_select = 1234))
x2 <- purrr::map_df(1:10, ~train_models(.x, sim_train = sim2_train, sim_test = sim2_test, sim_beta_0 = sim2_beta_0, coeff_names = coeff_names2, seed_select = 567))

# x1 %>%
#   group_by(name) %>% 
#   summarise(across(starts_with("res"), mean))
# x1 %>%
#   group_by(name) %>% 
#   summarise(across(starts_with("res"), sd))
x2 %>%
  group_by(name) %>%
  summarise(across(starts_with("res"), mean))
x2 %>%
  group_by(name) %>%
  summarise(across(starts_with("res"), sd))

################################################################################
################################################################################
# Testing Ground
################################################################################

i = 1
x <- purrr::map_df(1:100, function(i){
  print(i)

  train_x <- as.matrix(sim1_train[[i]][,coeff_names1])
  train_y <- as.matrix(sim1_train[[i]][,"y"])

  test_x_1 <- as.matrix(sim1_test[[i]][,c("p_0",coeff_names1)])
  test_x <- as.matrix(sim1_test[[i]][,coeff_names1])
  test_y <- as.matrix(sim1_test[[i]][,"y"])

  # return(eval_metrics)
})