################################################################################
rm(list = ls())

setwd("~/My Drive/GitHub/RobZS/R")
source("allprogs_RobZS.R")
setwd("~/Downloads/GitHub/stat548-proj3/R")

source("0a_functions_utils.R")
# source("0c_functions_models.R")

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

################################################################################
n <- 154
p <- ncol(selbal_HIV_x)
rm(selbal_HIV_x)
# Training on full model
train_x <- selbal_train %>%
  slice(1:154) %>%
  dplyr::select(-c(row_index, HIV_Status)) %>%
  as.matrix()

train_y <- selbal_train %>%
  slice(1:154) %>%
  dplyr::pull(HIV_Status)
################################################################################

model_lasso <- readRDS("results/glm_lasso.RDS")
model_enetLTS <- readRDS("results/glm_enetLTS.RDS")
model_zeroSum <- readRDS("results/glm_zeroSum.RDS")
model_RobZS <- readRDS("results/glm_RobZS.RDS")

Z <- cbind(1, train_x)
beta_zeroSum <- coefficients(model_zeroSum)
beta_RobZS <- coefficients(model_RobZS)

deviances_RobZS <- compute_d_z_beta_y(y = train_y, beta = beta_RobZS, z = Z)
deviances_zeroSum <- compute_d_z_beta_y(y = train_y, beta = beta_zeroSum, z = Z)

score_RobZS <- Z%*%beta_RobZS
score_zeroSum <- Z%*%beta_zeroSum
################################################################################

df_RobZS <- as_tibble(list(type = "RobZS", 
                      deviance = as.vector(deviances_RobZS),
                      score = as.vector(score_RobZS),
                      outlier = as.numeric(model_RobZS$raw.wt == 0),
                      RobZS_outlier = as.numeric(model_RobZS$raw.wt == 0),
                      y_true = train_y)) %>%
  mutate(row_id = row_number())
df_zeroSum <- as_tibble(list(type = "zeroSum", 
                      deviance = as.vector(deviances_zeroSum),
                      score = as.vector(score_zeroSum),
                      # outlier = as.numeric(model_RobZS$raw.wt == 0),
                      RobZS_outlier = as.numeric(model_RobZS$raw.wt == 0),
                      y_true = train_y)) %>%
  mutate(row_id = row_number())

df_full <- bind_rows(df_RobZS, df_zeroSum) %>%
  # mutate(outlier = as.factor(outlier))
  # mutate(outlier = ifelse(is.na(outlier), 2, outlier),
  mutate(outlier = case_when(outlier == 0 ~ "No Outlier",
                             outlier == 1 ~ "Outlier",
                             is.na(outlier) ~ "Other")) %>%
  mutate(outlier = forcats::fct_relevel(outlier, "No Outlier", "Outlier", "Other"),
         RobZS_outlier = factor(RobZS_outlier),
         y_true = as.factor(y_true))

################################################################################

p <- ggplot(data = df_full) +
  geom_point(aes(x = score, y = deviance, shape = type, col = outlier), alpha = 0.65, size = 2, stroke = 1) +
  # scale_color_discrete(palette = "Set2") +
  scale_color_brewer(palette = "Set1") +
  scale_shape_manual(values = c(2, 5)) +
  labs(x = "Score", y = "Deviance", col = "Type", shape = "Model") +
  theme_light() +
  theme(legend.position = "bottom")

ggsave(p, filename = "figures/selbal_HIV_analysis_deviance.pdf", width = 6, height = 5)

df_plot <- df_full %>% 
  dplyr::select(RobZS_outlier, row_id, type, score, y_true) %>%
  pivot_wider(names_from = c("type"), values_from = "score")

p <- ggplot(df_plot) +
  geom_point(aes(x = RobZS, y = zeroSum, shape = y_true, col = RobZS_outlier), alpha = 1, stroke = 1) +
  scale_color_brewer(palette = "Set2") +
  scale_shape_manual(values = c(1, 3)) +
  labs(shape = "HIV Positive", col = "RobZS Outlier") +
  theme_light() +
  theme(legend.position = "bottom")

ggsave(p, filename = "figures/selbal_HIV_analysis_outliers.pdf", width = 6, height = 5)

################################################################################
# Comparison with models after removing the outliers

# model_lasso_no_outliers <- readRDS("results/glm_lasso_no_outliers.RDS")
# model_enetLTS_no_outliers <- readRDS("results/glm_enetLTS_no_outliers.RDS")
# model_zeroSum_no_outliers <- readRDS("results/glm_zeroSum_no_outliers.RDS")
# model_RobZS_no_outliers <- readRDS("results/glm_RobZS_no_outliers.RDS")

################################################################################
compute_eval_new <- function(preds, train_y, name, threshold = 0.5){
  pred_hard <- as.numeric(preds > threshold)
  res_auc <- as.numeric(pROC::roc(response = as.vector(train_y), predictor = as.vector(preds), quiet = TRUE)$auc)
  res_se <- caret::sensitivity(as.factor(pred_hard), as.factor(train_y))
  res_sp <- caret::specificity(as.factor(pred_hard), as.factor(train_y))
  
  return(as_tibble(list(res_se = res_se, res_sp = res_sp, res_auc = res_auc)) %>%
           mutate(type = name))
}

outlier_id <- which(model_RobZS$raw.wt == 0)

train_y_new <- train_y[-outlier_id]

y_pred_LASSO <- predict(model_lasso, newx = train_x, type = "response")
y_pred_LASSO_new <- y_pred_LASSO[-outlier_id]

y_pred_enetLTS <- predict(model_enetLTS, newX = train_x, type = "response", vers = "reweighted")$reweighted.response
y_pred_enetLTS_new <- y_pred_enetLTS[-outlier_id]

y_pred_zeroSum <- compute_mu_beta_z(beta = coefficients(model_zeroSum), z = cbind(1,train_x))
y_pred_zeroSum_new <- y_pred_zeroSum[-outlier_id]

y_pred_RobZS <- model_RobZS$fitted.values
y_pred_RobZS_new <- model_RobZS$fitted.values[-outlier_id]


threshold <- 0.5

pred_compare_0.5 <- bind_rows(
  compute_eval_new(preds = y_pred_LASSO, train_y = train_y, name = "LASSO", threshold = threshold),
  compute_eval_new(preds = y_pred_LASSO_new, train_y = train_y_new, name = "LASSO_new", threshold = threshold),
  compute_eval_new(preds = y_pred_enetLTS, train_y = train_y, name = "RobLL", threshold = threshold),
  compute_eval_new(preds = y_pred_enetLTS_new, train_y = train_y_new, name = "RobLL_new", threshold = threshold),
  compute_eval_new(preds = y_pred_zeroSum, train_y = train_y, name = "zeroSum", threshold = threshold),
  compute_eval_new(preds = y_pred_zeroSum_new, train_y = train_y_new, name = "zeroSum_new", threshold = threshold),
  compute_eval_new(preds = y_pred_RobZS, train_y = train_y, name = "RobZS", threshold = threshold),
  compute_eval_new(preds = y_pred_RobZS_new, train_y = train_y_new, name = "RobZS_new", threshold = threshold)) %>%
  relocate(type)

readr::write_csv(pred_compare_0.5, "results/selbal_pred_compare_0_5.csv")

threshold <- 0.75

pred_compare_0.75 <- bind_rows(
  compute_eval_new(preds = y_pred_LASSO, train_y = train_y, name = "LASSO", threshold = threshold),
  compute_eval_new(preds = y_pred_LASSO_new, train_y = train_y_new, name = "LASSO_new", threshold = threshold),
  compute_eval_new(preds = y_pred_enetLTS, train_y = train_y, name = "RobLL", threshold = threshold),
  compute_eval_new(preds = y_pred_enetLTS_new, train_y = train_y_new, name = "RobLL_new", threshold = threshold),
  compute_eval_new(preds = y_pred_zeroSum, train_y = train_y, name = "zeroSum", threshold = threshold),
  compute_eval_new(preds = y_pred_zeroSum_new, train_y = train_y_new, name = "zeroSum_new", threshold = threshold),
  compute_eval_new(preds = pred_RobZS, train_y = train_y, name = "RobZS", threshold = threshold),
  compute_eval_new(preds = pred_RobZS_new, train_y = train_y_new, name = "RobZS_new", threshold = threshold)) %>%
  relocate(type)

readr::write_csv(pred_compare_0.75, "results/selbal_pred_compare_0_75.csv")

################################################################################

df <- bind_rows(
  as_tibble(list(y_pred = y_pred_LASSO, y_true = train_y, name = "LASSO", outlier = model_RobZS$raw.wt == 0)),
  as_tibble(list(y_pred = y_pred_enetLTS, y_true = train_y, name = "RobLL", outlier = model_RobZS$raw.wt == 0)),
  as_tibble(list(y_pred = y_pred_zeroSum, y_true = train_y, name = "ZeroSum", outlier = model_RobZS$raw.wt == 0)),
  as_tibble(list(y_pred = y_pred_RobZS, y_true = train_y, name = "RobZS", outlier = model_RobZS$raw.wt == 0))
  ) %>%
  mutate(outlier = as.factor(as.numeric(outlier)))

df_new <- bind_rows(
  as_tibble(list(y_pred = y_pred_LASSO_new, y_true = train_y_new, name = "LASSO")),
  as_tibble(list(y_pred = y_pred_enetLTS_new, y_true = train_y_new, name = "RobLL")),
  as_tibble(list(y_pred = y_pred_zeroSum_new, y_true = train_y_new, name = "ZeroSum")),
  as_tibble(list(y_pred = y_pred_RobZS_new, y_true = train_y_new, name = "RobZS"))
  )

# df_full <-bind_rows(
#   as_tibble(list(y_pred = y_pred_LASSO, y_true = train_y, name = "LASSO", outlier = model_RobZS$raw.wt == 0)),
#   as_tibble(list(y_pred = y_pred_enetLTS, y_true = train_y, name = "RobLL", outlier = model_RobZS$raw.wt == 0)),
#   as_tibble(list(y_pred = y_pred_zeroSum, y_true = train_y, name = "ZeroSum", outlier = model_RobZS$raw.wt == 0)),
#   as_tibble(list(y_pred = y_pred_RobZS, y_true = train_y, name = "RobZS", outlier = model_RobZS$raw.wt == 0)),
#   as_tibble(list(y_pred = y_pred_LASSO_new, y_true = train_y_new, name = "LASSO (new)", outlier = 0)),
#   as_tibble(list(y_pred = y_pred_enetLTS_new, y_true = train_y_new, name = "RobLL (new)", outlier = 0)),
#   as_tibble(list(y_pred = y_pred_zeroSum_new, y_true = train_y_new, name = "ZeroSum (new)", outlier = 0)),
#   as_tibble(list(y_pred = y_pred_RobZS_new, y_true = train_y_new, name = "RobZS (new)", outlier = 0))
# ) %>%
#   mutate(outlier = as.factor(as.numeric(outlier)))

p1 <- ggplot(df) +
  # geom_point(aes(x = y_pred, y = y_true, col = outlier), position = "jitter") +
  geom_jitter(aes(x = y_pred, y = y_true, col = outlier), width = 0, height = 0.1, alpha = 0.5) +
  facet_wrap(~name, nrow = 2) +
  expand_limits(x = c(0,1)) +
  labs(x = "Fitted Response", y = "True Response", col = "Outlier (RobZS)") +
  theme_light() +
  theme(legend.position = "bottom")

ggsave(p1, filename = "figures/fitted_values.pdf", width = 6, height = 5)

p2 <- ggplot(df_new) +
  geom_jitter(aes(x = y_pred, y = y_true), col = "red", width = 0, height = 0.1, alpha = 0.5) +
  facet_wrap(~name) +
  expand_limits(x = c(0,1)) +
  labs(x = "Fitted Response", y = "True Response") +
  theme_light()

ggsave(p2, filename = "figures/fitted_values_no_outliers.pdf", width = 6, height = 5)

