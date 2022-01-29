setwd("~/My Drive/GitHub/RobZS/R")
setwd("~/Downloads/GitHub/stat548-proj3/R")
rm(list = ls())

source("0b_functions_utils.R")


# zeroSum
# enetLTS
# glmnet

################################################################################

library(tidyverse)
library(glmnet)

# sim_all_1 <- readRDS("data/1_sim_all.RDS")
# Z_contaminated <- sim_all_1$Z_contaminated

sim1 <- readRDS("data/sim2.RDS")

sim1_train <- do.call("rbind", lapply(sim1, '[', "Z_contaminated"))
sim1_test <- do.call("rbind", lapply(sim1, '[', "Z_test"))
p <- length(sim1[[1]]$beta) - 1
n <- nrow(sim1[[1]]$Z_contaminated)

################################################################################

coeff_names <- purrr::map_chr(1:p, ~paste0("p_",.x))

res_mse = vector(mode = "numeric")
res_mae = vector(mode = "numeric")
res_ml = vector(mode = "numeric")
res_se = vector(mode = "numeric")
res_sp = vector(mode = "numeric")
res_auc = vector(mode = "numeric")

################################################################################

for(i in (1:100)){
  print(i)
  test_df_x <- sim1_train[[i]][,coeff_names] %>%
    as.matrix()
  
  test_df_y <- sim1_train[[i]][,"y"] %>%
    as.matrix()
  
  test_df_x_new <- sim1_test[[i]][,coeff_names] %>%
    as.matrix()
  test_df_x_new_1 <- sim1_test[[i]][,c("p_0",coeff_names)] %>%
    as.matrix()
  test_df_x_new_y <- sim1_test[[i]][,"y"] %>%
    as.matrix()
  
  glm1 <- glmnet(x = test_df_x, y = test_df_y,
                 family = "binomial", 
                 lambda = cv.glmnet(x = test_df_x, y = test_df_y)$lambda.1se,
                 alpha = 1)
  
  y1 <- glmnet::predict.glmnet(glm1, newx = test_df_x_new)
  
  new_error <- compute_mu_beta_z(beta = coefficients(glm1), z = test_df_x_new_1) - test_df_x_new_y
  res_mse <- c(res_mse, mean(new_error^2))
  res_mae <- c(res_mae, mean(abs(new_error)))
}

min(res_mse)
min(res_mae)
sd(res_mse)
sd(res_mae)

# X1 <- Z_contaminated[,coeff_names] %>%
#   as.matrix()
# 
# Y1 <- Z_contaminated[,"y"] %>%
#   as.matrix()
# 
# test_robzs <- RobZS(xx = X1, yy = Y1, family = "binomial")
# coef.RobZS(test_robzs)
