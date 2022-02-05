################################################################################
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
library(cvTools)

dir.create("results")
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

x1 <- purrr::map_df(1:100, ~train_models(.x, sim_train = sim1_train, sim_test = sim1_test, sim_beta_0 = sim1_beta_0, coeff_names = coeff_names1, seed_select = 223))
# readr::write_csv(x1, "results/sim1_results.csv")
x1 <- read_csv("results/sim1_results.csv")

x2 <- purrr::map_df(1:100, ~train_models(.x, sim_train = sim2_train, sim_test = sim2_test, sim_beta_0 = sim2_beta_0, coeff_names = coeff_names2, seed_select = 567))
# readr::write_csv(x2, "results/sim2_results.csv")
x2 <- read_csv("results/sim2_results.csv")

x1 %>%
  group_by(name) %>%
  summarise(across(starts_with("res"), list(mean = mean,
                                            sd = sd), na.rm = TRUE)) %>%
  relocate(name, res_se_mean, res_se_sd, res_sp_mean, res_sp_sd, res_auc_mean, res_auc_sd) %>%
  readr::write_csv("results/table_sim1.csv")
x2 %>%
  group_by(name) %>%
  summarise(across(starts_with("res"), list(mean = mean,
                                            sd = sd), na.rm = TRUE)) %>%
  relocate(name, res_se_mean, res_se_sd, res_sp_mean, res_sp_sd, res_auc_mean, res_auc_sd) %>%
  readr::write_csv("results/table_sim2.csv")
# x1 %>%
#   group_by(name) %>%
#   summarise(across(starts_with("res"), sd, na.rm = TRUE))
# x2 %>%
#   group_by(name) %>%
#   summarise(across(starts_with("res"), mean))
# x2 %>%
#   group_by(name) %>%
#   summarise(across(starts_with("res"), sd))

################################################################################
################################################################################
# Testing Ground
################################################################################
# 
# i = 1
# x <- purrr::map_df(1:100, function(i){
#   print(i)
# 
#   train_x <- as.matrix(sim1_train[[i]][,coeff_names1])
#   train_y <- as.matrix(sim1_train[[i]][,"y"])
# 
#   test_x_1 <- as.matrix(sim1_test[[i]][,c("p_0",coeff_names1)])
#   test_x <- as.matrix(sim1_test[[i]][,coeff_names1])
#   test_y <- as.matrix(sim1_test[[i]][,"y"])
#   
#   sim_beta_0 <- sim1_beta_0
#   seed = 5
# 
#   # return(eval_metrics)
# })