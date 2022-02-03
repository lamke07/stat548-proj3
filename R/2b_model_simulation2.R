rm(list = ls())

setwd("~/My Drive/GitHub/RobZS/R")
source("allprogs_RobZS.R")
setwd("~/Downloads/GitHub/stat548-proj3/R")

source("0a_functions_utils.R")
source("0c_functions_models.R")

library(tidyverse)
library(glmnet)
library(enetLTS) # RobLL
library(zeroSum) # LZS
library(caret)
library(pROC)

################################################################################
sim1 <- readRDS("data/extended_sim1.RDS")

sim1_train <- do.call("rbind", lapply(sim1, '[', "Z_contaminated"))
sim1_beta <- do.call("rbind", lapply(sim1, '[', "beta"))[[1]]
sim1_test <- do.call("rbind", lapply(sim1, '[', "Z_test"))

sim1_beta_0 <- sim1_beta[-1]
coeff_names1 <- purrr::map_chr(1:length(sim1_beta_0), ~paste0("p_",.x))

################################################################################

sim2 <- readRDS("data/extended_sim2.RDS")

sim2_train <- do.call("rbind", lapply(sim2, '[', "Z_contaminated"))
sim2_beta <- do.call("rbind", lapply(sim2, '[', "beta"))[[1]]
sim2_test <- do.call("rbind", lapply(sim2, '[', "Z_test"))

sim2_beta_0 <- sim2_beta[-1]
coeff_names2 <- purrr::map_chr(1:length(sim2_beta_0), ~paste0("p_",.x))


################################################################################
################################################################################

train_models <- function(i, sim_train, sim_test, sim_beta_0, coeff_names, seed_select = 123, ncores = 1){
  print(i)
  # Preprocessing
  train_x <- as.matrix(sim_train[[i]][,coeff_names])
  train_y <- as.matrix(sim_train[[i]][,"y"])
  
  test_x_1 <- as.matrix(sim_test[[i]][,c("p_0",coeff_names)])
  test_x <- as.matrix(sim_test[[i]][,coeff_names])
  test_y <- as.matrix(sim_test[[i]][,"y"])
  
  seed <- i + seed_select + length(sim_beta_0) + nrow(train_x)
  
  # Lasso estimates
  eval_metrics_lasso <- fit_lasso(train_x = train_x, train_y = train_y, 
                                  test_x = test_x, test_y = test_y,
                                  seed = seed, sim_beta_0 = sim_beta_0)
  
  # Lasso estimates
  eval_metrics_LTS <- fit_enetLTS(train_x = train_x, train_y = train_y, 
                                  test_x = test_x, test_y = test_y,
                                  seed = seed, sim_beta_0 = sim_beta_0, ncores = 6)
  
  # LZS estimates
  eval_metrics_LZS <- fit_zeroSum(train_x = train_x, train_y = train_y, 
                                  test_x = test_x_1, test_y = test_y,
                                  seed = seed, sim_beta_0 = sim_beta_0)
  
  eval_metrics_RobZS <- fit_RobZS(train_x = train_x, train_y = train_y, 
                                  test_x = test_x, test_y = test_y,
                                  seed = seed, sim_beta_0 = sim_beta_0)
  
  
  return(bind_rows(eval_metrics_lasso, eval_metrics_LTS, eval_metrics_LZS, eval_metrics_RobZS))
}

################################################################################


x1 <- purrr::map_df(1:2, ~train_models(.x, sim_train = sim1_train, sim_test = sim1_test, sim_beta_0 = sim1_beta_0, coeff_names = coeff_names1, seed_select = 1234))
# x2 <- purrr::map_df(1:10, ~train_models(.x, sim_train = sim2_train, sim_test = sim2_test, sim_beta_0 = sim2_beta_0, coeff_names = coeff_names2, seed_select = 567))

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