setwd("~/My Drive/GitHub/RobZS/R")
setwd("~/Downloads/GitHub/stat548-proj3/R")

rm(list = ls())

# zeroSum
# enetLTS
# glmnet

################################################################################

library(tidyverse)
library(glmnet)

# sim_all_1 <- readRDS("data/1_sim_all.RDS")
# Z_contaminated <- sim_all_1$Z_contaminated

sim1 <- readRDS("data/sim1.RDS")

sim1_data <- do.call("rbind", lapply(sim1, '[', "Z_contaminated"))
sim1_beta <- sim1[[1]]$beta 

################################################################################

p <- length(sim1_beta) - 1
coeff_names <- purrr::map_chr(1:p, ~paste0("p_",.x))

test_df_x <- sim1_data[[1]][,coeff_names] %>%
  as.matrix()

test_df_y <- sim1_data[[1]][,"y"] %>%
  as.matrix()

glm1 <- glmnet(x = test_df, y = test_df_y, family = "binomial", lambda = cv.glmnet(x = test_df, y = test_df_y)$lambda.1se, alpha = 1)

# X1 <- Z_contaminated[,coeff_names] %>%
#   as.matrix()
# 
# Y1 <- Z_contaminated[,"y"] %>%
#   as.matrix()
# 
# test_robzs <- RobZS(xx = X1, yy = Y1, family = "binomial")
# coef.RobZS(test_robzs)
