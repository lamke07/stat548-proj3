# setwd("~/My Drive/GitHub/RobZS/R")
# setwd("~/Downloads/GitHub/stat548-proj3/R")
rm(list = ls())

source("0b_functions_utils.R")

library(tidyverse)
library(glmnet)
library(caret)
library(pROC)
library(enetLTS)

# RobZs
# zeroSum
# enetLTS
# glmnet

################################################################################
sim1 <- readRDS("data/sim1.RDS")

sim1_train <- do.call("rbind", lapply(sim1, '[', "Z_contaminated"))
sim1_beta <- do.call("rbind", lapply(sim1, '[', "beta"))[[1]]
sim1_test <- do.call("rbind", lapply(sim1, '[', "Z_test"))

sim2 <- readRDS("data/sim2.RDS")

sim2_train <- do.call("rbind", lapply(sim2, '[', "Z_contaminated"))
sim2_beta <- do.call("rbind", lapply(sim2, '[', "beta"))[[1]]
sim2_test <- do.call("rbind", lapply(sim2, '[', "Z_test"))


################################################################################

sim1_beta_0 <- sim1_beta[-1]
coeff_names1 <- purrr::map_chr(1:length(sim1_beta_0), ~paste0("p_",.x))
sim2_beta_0 <- sim2_beta[-1]
coeff_names2 <- purrr::map_chr(1:length(sim2_beta_0), ~paste0("p_",.x))

################################################################################

train_models <- function(i, sim_train, sim_test, sim_beta_0, coeff_names, seed_select = 123){
  print(i)
  
  train_x <- as.matrix(sim_train[[i]][,coeff_names])
  train_y <- as.matrix(sim_train[[i]][,"y"])
  
  test_x <- as.matrix(sim_test[[i]][,coeff_names])
  test_x_1 <- as.matrix(sim_test[[i]][,c("p_0",coeff_names)])
  test_y <- as.matrix(sim_test[[i]][,"y"])
  
  set.seed(123 + i + seed_select + length(sim1_beta_0) + nrow(train_x))
  lambda_lasso <- cv.glmnet(x = train_x, y = train_y)$lambda.1se
  glm_lasso <- glmnet(x = train_x, y = train_y,
                      family = "binomial", 
                      lambda = lambda_lasso,
                      alpha = 1)
  
  y_pred <- predict(glm_lasso, newx = test_x, type = "response")
  
  eval_metrics <- compute_metrics(y_pred = y_pred, beta_pred = glm_lasso$beta,
                                  test_y = test_y, true_beta = sim_beta_0,
                                  name = "LASSO")
  
  return(eval_metrics)
}

x1 <- purrr::map_df(1:100, ~train_models(.x, sim_train = sim1_train, sim_test = sim1_test, sim_beta_0 = sim1_beta_0, coeff_names = coeff_names1, seed_select = 1234))
x2 <- purrr::map_df(1:100, ~train_models(.x, sim_train = sim2_train, sim_test = sim2_test, sim_beta_0 = sim2_beta_0, coeff_names = coeff_names2, seed_select = 567))

x1 %>%
  group_by(name) %>% 
  summarise(across(starts_with("res"), mean))
x1 %>%
  group_by(name) %>% 
  summarise(across(starts_with("res"), sd))
x2 %>%
  group_by(name) %>% 
  summarise(across(starts_with("res"), mean))
x2 %>%
  group_by(name) %>% 
  summarise(across(starts_with("res"), sd))


i = 1
x <- purrr::map_df(1:100, function(i){
  print(i)

  train_x <- as.matrix(sim1_train[[i]][,coeff_names1])
  train_y <- as.matrix(sim1_train[[i]][,"y"])

  test_x <- as.matrix(sim1_test[[i]][,coeff_names1])
  test_x_1 <- as.matrix(sim1_test[[i]][,c("p_0",coeff_names1)])
  test_y <- as.matrix(sim1_test[[i]][,"y"])

  # set.seed(123 + i + seed_select + length(sim1_beta_0) + nrow(train_x))
  # lambda_lasso <- cv.glmnet(x = train_x, y = train_y)$lambda.1se
  # glm_lasso <- glmnet(x = train_x, y = train_y,
  #                family = "binomial",
  #                lambda = lambda_lasso,
  #                alpha = 1)
  # 
  # y_pred <- predict(glm_lasso, newx = test_x, type = "response")
  
  glm_enetLTS <- enetLTS(xx = train_x, yy = as.vector(train_y), family = "binomial", alphas = 1, ncores = 6, plot = "FALSE")
  y_pred <- predict(a, newX = test_x, type = "response", vers = "reweighted")$reweighted.response

  eval_metrics <- compute_metrics(y_pred = y_pred, beta_pred = glm_enetLTS$coefficients,
                                  test_y = test_y, true_beta = sim1_beta_0,
                                  name = "LASSO")

  return(eval_metrics)
})

# a <- enetLTS(xx = train_x, yy = as.vector(train_y), family = "binomial", alphas = 1, ncores = 6, plot = "FALSE")
# b <- predict(a, newX = test_x, type = "response", vers = "reweighted")

x %>%
  group_by(name) %>% 
  summarise(across(starts_with("res"), mean))
x %>%
  group_by(name) %>% 
  summarise(across(starts_with("res"), sd))
