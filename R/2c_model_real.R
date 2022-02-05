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
                            dplyr::select(-row_sum) %>% 
                            log(),
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

# selbal_HIV %>%
#   dplyr::select(-MSM) %>%
#   group_by(HIV_Status) %>%
#   summarise(across(everything(), max)) %>%
#   ungroup() %>%
#   pivot_longer(-HIV_Status, names_to = "type", values_to = "max_count") %>%
#   # mutate(count_group = cut(max_count, breaks = c(0, 50, 100, 200, 300, 600, 800, 1000, 1500, 6000), labels = c("< 0050", "< 0100", "< 0200", "< 0300", "< 0600", "< 0800", "< 1000", "< 1500", "< 6000"))) %>%
#   mutate(count_group = cut(max_count, breaks = c(0, 100, 300, 600, 1000, 6000), labels = c("< 0100", "< 0300", "< 0600", "< 1500", "< 6000"))) %>%
#   group_by(HIV_Status, count_group) %>%
#   summarise(count = n())

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
                  cat("\n Replication:", fold_select, "| Fold: ", i, "\n")
                  
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
                  
                  res_new <- train_models_CV(train_x = train_x, train_y = train_y, 
                                             test_x = test_x, test_y = test_y, standardize = TRUE,
                                             seed = i*fold_select*100 + 5, ncores = 6) %>%
                    mutate(fold = fold_select, fold_test = i)
                  
                  return(res_new)
                })

  
  return(res)
}

train_models_CV <- function(train_x, train_y, test_x, test_y, sim_beta_0 = FALSE, seed = 123, ncores = ncores, standardize = TRUE){
  
  print("Fitting Lasso...")
  eval_metrics_lasso <- fit_lasso(train_x = train_x, train_y = train_y, 
                                  test_x = test_x, test_y = test_y,
                                  seed = seed, sim_beta_0 = sim_beta_0, standardize = standardize)
  
  print("Fitting LTS (RobLL)...")
  eval_metrics_LTS <- fit_enetLTS(train_x = train_x, train_y = train_y,
                                  test_x = test_x, test_y = test_y,
                                  seed = seed, sim_beta_0 = sim_beta_0, standardize = standardize,
                                  ncores = ncores)
  
  print("Fitting LZS...")
  eval_metrics_LZS <- fit_zeroSum(train_x = train_x, train_y = train_y, 
                                  test_x = cbind(1, test_x), test_y = test_y, standardize = standardize,
                                  seed = seed, sim_beta_0 = sim_beta_0)
  
  
  # print("Fitting RobZS...")
  # eval_metrics_RobZS <- fit_RobZS(train_x = train_x, train_y = train_y,
  #                                 test_x = test_x, test_y = test_y,
  #                                 seed = seed, sim_beta_0 = sim_beta_0, standardize = standardize,
  #                                 ncores = ncores)
  
  # return(bind_rows(eval_metrics_lasso, eval_metrics_LTS, eval_metrics_LZS, eval_metrics_RobZS))
  return(bind_rows(eval_metrics_lasso, eval_metrics_LTS, eval_metrics_LZS))
  # return(bind_rows(eval_metrics_lasso, eval_metrics_LZS))
}

################################################################################
safe_train_models_all <- safely(.f = train_models_all)
# x <- train_models_all(df = selbal_train, fold_select = 1, cv_folds = cv_folds)
y <- purrr::map(1:20, ~safe_train_models_all(df = selbal_train, fold_select = .x, cv_folds = cv_folds))
y_result <- purrr::map_df(x1, "result")
y_error <- purrr::map(x1, "error")
readr::write_csv(y_result, "results/sim1_results.csv")
# saveRDS(y, "results/selbal_y.RDS")

dir.create("figures")

# Average coefficient

df <- y %>% 
  group_by(name) %>% 
  summarise(across(starts_with("p_"), mean)) %>%
  ungroup() %>% 
  pivot_longer(cols= starts_with("p_"), names_to = "beta_coef", values_to = "beta_val") 

p <- ggplot(data = df) +
  geom_point(aes(x = beta_coef, y = beta_val, col = name)) +
  labs(x = "Component", y = "Coefficient for component", col = "Model", title = "Average coefficient over 100 CV folds") +
  theme_light() +
  theme(axis.text.x = element_text(angle = 90))

ggsave(p, filename = "figures/selbal_HIV_average_coefficient.pdf", width = 10, height = 6)

# Evaluation metrics

df <- y %>%
  dplyr::select(name, starts_with("res")) %>%
  pivot_longer(cols = starts_with("res"), names_to = "res_type", values_to = "res")

p <- ggplot(data = df) +
  geom_boxplot(aes(x = name, y = res, group = name)) +
  facet_wrap(~res_type) +
  labs(x = "Models", y = "Evaluation Metric", title = "Comparison of Evaluation metrics for models") +
  theme_light()

ggsave(p, filename = "figures/selbal_HIV_average_evaluation.pdf", width = 10, height = 6)

# Proportion of zeros

df <- y %>%
  count(name, no_zeros) %>%
  arrange(name, -no_zeros) %>%
  group_by(name) %>%
  mutate(cumsum_no_zeros = cumsum(n),
         total_zeros = sum(n)) %>%
  ungroup() %>%
  arrange(name, no_zeros) %>%
  mutate(frac_no_zeros = cumsum_no_zeros/total_zeros)

p <- ggplot(df) +
  geom_line(aes(x = no_zeros, y = frac_no_zeros, colour = factor(name))) +
  labs(y = "P(|{beta == 0}| > #zeros)", x = "#zeros", colour = "Model",
       title = "Proportion of coefficients with at least a given number of zeros") +
  theme_light()

ggsave(p, filename = "figures/selbal_HIV_proportion.pdf", width = 6, height = 4)

# ################################################################################
# ################################################################################
# ################################################################################
# train_x <- selbal_train %>%
#   dplyr::select(-c(row_index, HIV_Status)) %>%
#   as.matrix()
# 
# train_y <- selbal_train %>%
#   dplyr::pull(HIV_Status)
# 
# ################################################################################
# set.seed(123)
# lambda_lasso <- cv.glmnet(x = train_x, y = train_y, alpha = 1, family = "binomial")$lambda.min
# glm_lasso <- glmnet(x = train_x, y = train_y, family = "binomial", lambda = lambda_lasso, alpha = 1)
# saveRDS(glm_lasso, "results/glm_lasso.RDS")
# 
# # # Obtain coefficients
# # glm_lasso_beta <- glm_lasso$beta
# # beta_full <- t(as.vector(coef(glm_lasso)))
# # colnames(beta_full) = purrr::map_chr(0:ncol(train_x), ~paste0("p_",.x))
# 
# # y_pred_lasso <- predict(glm_lasso, newx = test_x, type = "response")
# ################################################################################
# glm_enetLTS <- enetLTS(xx = train_x, yy = as.vector(train_y), family = "binomial", alphas = 1, ncores = 6, seed = 234, plot = "FALSE")
# saveRDS(glm_enetLTS, "results/glm_enetLTS.RDS")
# 
# # # Obtain coefficients
# # glm_enetLTS_beta <- glm_enetLTS$coefficients
# # beta_full <- t(as.vector(coefficients(glm_enetLTS)))
# # colnames(beta_full) = purrr::map_chr(0:ncol(train_x), ~paste0("p_",.x))
# 
# # y_pred_LTS <- predict(glm_enetLTS, newX = test_x, type = "response", vers = "reweighted")$reweighted.response
# ################################################################################
# # LZS estimates
# set.seed(245)
# glm_zeroSum <- zeroSum(train_x, as.vector(train_y), family = "binomial", alpha = 1)
# saveRDS(glm_zeroSum, "results/glm_zeroSum.RDS")
# 
# # # Obtain coefficients
# # glm_zeroSum_beta <- coefficients(glm_zeroSum, s = "lambda.min")
# # beta_full <- t(as.vector(glm_zeroSum_beta))
# # colnames(beta_full) = purrr::map_chr(0:ncol(train_x), ~paste0("p_",.x))
# 
# # # Obtain predictions
# # y_pred_zeroSum <- compute_mu_beta_z(beta = glm_zeroSum_beta, z = test_x)
# ################################################################################
# glm_RobZS <- RobZS(xx = train_x, yy = as.vector(train_y), family = "binomial", alphas = 1, seed = 345, ncores = 6, plot = FALSE)
# saveRDS(glm_RobZS, "results/glm_RobZS.RDS")
# 
# # # Obtain coefficients
# # glm_RobZS_beta <- glm_RobZS$coefficients
# # beta_full <- t(as.vector(coefficients(glm_RobZS)))
# # colnames(beta_full) = purrr::map_chr(0:ncol(train_x), ~paste0("p_",.x))
# 
# # # Obtain predictions
# # y_pred_RobZS <- predict.RobZS(glm_RobZS, newX = test_x, vers = "reweighted", type0 = "response", intercept = TRUE)$reweighted.response
# ################################################################################
# 
# 
# 
