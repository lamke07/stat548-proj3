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
sim1 <- readRDS("data/extended2_sim1.RDS")

sim1_train <- do.call("rbind", lapply(sim1, '[', "Z_contaminated"))
sim1_beta <- do.call("rbind", lapply(sim1, '[', "beta"))[[1]]
sim1_test <- do.call("rbind", lapply(sim1, '[', "Z_test"))

sim1_beta_0 <- sim1_beta[-1]
coeff_names1 <- purrr::map_chr(1:length(sim1_beta_0), ~paste0("p_",.x))

################################################################################

sim2 <- readRDS("data/extended2_sim2.RDS")

sim2_train <- do.call("rbind", lapply(sim2, '[', "Z_contaminated"))
sim2_beta <- do.call("rbind", lapply(sim2, '[', "beta"))[[1]]
sim2_test <- do.call("rbind", lapply(sim2, '[', "Z_test"))

sim2_beta_0 <- sim2_beta[-1]
coeff_names2 <- purrr::map_chr(1:length(sim2_beta_0), ~paste0("p_",.x))

################################################################################
################################################################################

x1 <- purrr::map_df(1:100, ~train_models(.x, sim_train = sim1_train, sim_test = sim1_test, sim_beta_0 = sim1_beta_0, coeff_names = coeff_names1, seed_select = 3390))
# readr::write_csv(x1, "results/extended_sim1_results.csv")

x2 <- purrr::map_df(1:100, ~train_models(.x, sim_train = sim2_train, sim_test = sim2_test, sim_beta_0 = sim2_beta_0, coeff_names = coeff_names2, seed_select = 567))
readr::write_csv(x2, "results/extended_sim2_results.csv")

x1 %>%
  group_by(name) %>%
  summarise(across(starts_with("res"), list(mean = mean,
                                            sd = sd), na.rm = TRUE)) %>%
  relocate(name, res_se_mean, res_se_sd, res_sp_mean, res_sp_sd, res_auc_mean, res_auc_sd) %>%
  readr::write_csv("results/table_extended_sim1.csv")
