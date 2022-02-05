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

selbal_HIV <- readr::read_csv("data/selbal_HIV.csv") %>%
  mutate(HIV_Status = ifelse(HIV_Status == "Pos", 1, 0))

selbal_HIV_x <- selbal_HIV %>%
  dplyr::select(-starts_with("HIV_Status"), -starts_with("MSM")) %>%
  mutate(across(everything(), ~ifelse(.x == 0, 0.5, .x)))
# selbal_HIV_MSM <- data.frame(MSM = selbal_HIV[,61])

# Get training data
selbal_train <- bind_cols(normalize_rows(selbal_HIV_x) %>%
                            dplyr::select(-row_sum),
                          selbal_HIV %>%
                            dplyr::select(HIV_Status)) %>%
  arrange(HIV_Status) %>%
  mutate(row_index = row_number()) %>%
  relocate(row_index, HIV_Status)

n <- nrow(selbal_train)
n0 <- nrow(selbal_train %>% filter(HIV_Status == 0))
n1 <- nrow(selbal_train %>% filter(HIV_Status == 1))
p <- ncol(selbal_HIV_x)

rm(selbal_HIV_x)

# # No column with more than 90% zeros
# which(colSums(selbal_HIV == 0) > 135)
# # 128 HIV positive
# sum(selbal_HIV$HIV_Status)
# # 35% 0's
# sum(selbal_HIV %>% dplyr::select(-HIV_Status) == 0)/(155*60)

################################################################################
################################################################################
################################################################################
train_x <- selbal_train %>%
  dplyr::select(-c(row_index, HIV_Status)) %>%
  as.matrix()

train_y <- selbal_train %>%
  dplyr::pull(HIV_Status)

################################################################################
# Model runs
set.seed(123)
lambda_lasso <- cv.glmnet(x = train_x, y = train_y, alpha = 1, family = "binomial")$lambda.min
glm_lasso <- glmnet(x = train_x, y = train_y, family = "binomial", lambda = lambda_lasso, alpha = 1)
# saveRDS(glm_lasso, "results/glm_lasso.RDS")

glm_enetLTS <- enetLTS(xx = train_x, yy = as.vector(train_y), family = "binomial", alphas = 1, ncores = 6, seed = 234, plot = "FALSE")
# saveRDS(glm_enetLTS, "results/glm_enetLTS.RDS")

set.seed(245)
glm_zeroSum <- zeroSum(train_x, as.vector(train_y), family = "binomial", alpha = 1, ncores = 6)
# saveRDS(glm_zeroSum, "results/glm_zeroSum.RDS")

glm_RobZS <- RobZS(xx = train_x, yy = as.vector(train_y), family = "binomial", alphas = 1, seed = 345, ncores = 6, plot = FALSE)
# saveRDS(glm_RobZS, "results/glm_RobZS.RDS")

################################################################################
