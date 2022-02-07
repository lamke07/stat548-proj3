################################################################################
rm(list = ls())

setwd("~/My Drive/GitHub/RobZS/R")
source("allprogs_RobZS.R")
setwd("~/Downloads/GitHub/stat548-proj3/R")

source("0a_functions_utils.R")
source("0c_functions_models.R")

dir.create("results")

################################################################################

selbal_HIV <- readr::read_csv("data/selbal_HIV.csv") %>%
  mutate(HIV_Status = ifelse(HIV_Status == "Pos", 1, 0))

selbal_HIV_x <- selbal_HIV %>%
  dplyr::select(-starts_with("HIV_Status"), -starts_with("MSM")) %>%
  mutate(across(everything(), ~ifelse(.x == 0, 0.5, .x)))
# selbal_HIV_MSM <- data.frame(MSM = selbal_HIV[,61])

# Get training data
selbal_train <- bind_cols(selbal_HIV_x %>%
                            normalize_rows() %>%
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

################################################################################
################################################################################
# df <- selbal_train
# fold_select <- 1
# cv_folds <- cv_folds

train_models_all <- function(df, cv_folds, fold_select){
  cv_fold <- cv_folds[[fold_select]]
  n <- nrow(df)

  res <- purrr::map_df(1:3,
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

################################################################################
# Run all cross validation
cv_folds <- lapply(1:20, balanced_folds, n0 = n0, n1 = n1, K = 3)
safe_train_models_all <- safely(.f = train_models_all)

x_selbal <- purrr::map(1:20, ~safe_train_models_all(df = selbal_train, fold_select = .x, cv_folds = cv_folds))
x_selbal_result <- purrr::map_df(x_selbal, "result")
x_selbal_error <- purrr::map(x_selbal, "error")

saveRDS(x_selbal, "results/selbal_hiv_results.RDS")
readr::write_csv(x_selbal_result, "results/selbal_hiv_results.csv")
x_selbal_result %>%
  group_by(name) %>%
  summarise(across(starts_with("res"), list(mean = mean,
                                            sd = sd), na.rm = TRUE)) %>%
  relocate(name, res_se_mean, res_se_sd, res_sp_mean, res_sp_sd, res_auc_mean, res_auc_sd) %>%
  readr::write_csv("results/selbal_hiv_results-table.csv")

################################################################################

dir.create("figures")

# Average coefficient
df <- x_selbal_result %>% 
  dplyr::select(-p_0) %>%
  group_by(name) %>% 
  summarise(across(starts_with("p_"), mean)) %>%
  ungroup() %>% 
  pivot_longer(cols= starts_with("p_"), names_to = "beta_coef", values_to = "beta_val") 

p <- ggplot(data = df %>%
              mutate(beta_coef = as.numeric(substr(beta_coef, 3,5)))) +
  geom_point(aes(x = beta_coef, y = beta_val, col = name)) +
  labs(x = "Model Component", y = "Coefficient for component", col = "Model", title = "Average coefficient over 100 CV folds") +
  theme_light()

ggsave(p, filename = "figures/selbal_HIV_average_coefficient.pdf", width = 10, height = 6)

# Evaluation metrics
df <- x_selbal_result %>%
  dplyr::select(name, starts_with("res")) %>%
  pivot_longer(cols = starts_with("res"), names_to = "res_type", values_to = "res") %>%
  mutate(res_type = toupper(substr(res_type, 5, 9)))

p <- ggplot(data = df) +
  geom_boxplot(aes(x = name, y = res, group = name)) +
  facet_wrap(~res_type, scale = "free_y") +
  labs(x = "Models", y = "Evaluation Metric", title = "Comparison of Evaluation metrics for models") +
  theme_light()

ggsave(p, filename = "figures/selbal_HIV_average_evaluation.pdf", width = 10, height = 6)

# Proportion of zeros
df <- x_selbal_result %>%
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

################################################################################
################################################################################
################################################################################
# Training on full model
train_x <- selbal_train %>%
  slice(1:154) %>%
  dplyr::select(-c(row_index, HIV_Status)) %>%
  as.matrix()

train_y <- selbal_train %>%
  slice(1:154) %>%
  dplyr::pull(HIV_Status)

################################################################################
# Model runs
standardize = TRUE

set.seed(123)
lambda_lasso <- cv.glmnet(x = train_x, y = train_y, alpha = 1, family = "binomial", standardize = standardize)$lambda.min
glm_lasso <- glmnet(x = train_x, y = train_y, family = "binomial", lambda = lambda_lasso, alpha = 1, standardize = standardize)
# saveRDS(glm_lasso, "results/glm_lasso.RDS")

glm_enetLTS <- enetLTS(xx = train_x, yy = as.vector(train_y), family = "binomial", alphas = 1, ncores = 6, seed = 234, plot = "FALSE", scal = standardize)
# saveRDS(glm_enetLTS, "results/glm_enetLTS.RDS")

set.seed(245)
glm_zeroSum <- zeroSum(train_x, as.vector(train_y), family = "binomial", alpha = 1, ncores = 6, standardize = standardize)
# saveRDS(glm_zeroSum, "results/glm_zeroSum.RDS")

glm_RobZS <- RobZS(xx = train_x, yy = as.vector(train_y), family = "binomial", alphas = 1, seed = 567, ncores = 6, plot = FALSE, scal = standardize)
# saveRDS(glm_RobZS, "results/glm_RobZS.RDS")


