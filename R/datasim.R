################################################################################
rm(list = ls())

library(tidyverse)
library(MASS)
library(boot)

compute_mu_beta_z <- function(beta, z){
  # Input: vector[p] beta, matrix[n,p] z
  # Check all dimensions are the same
  stopifnot(dim(z)[2] == length(beta))
  inv.logit(as.numeric(z %*% beta))
}

################################################################################
# Hyperparameters

n <- 20
p <- 17
rho <- 0.2

coeff_names <- purrr::map_chr(0:p, ~paste0("p_",.x))

# Robust estimators
xi <- 3/4

# Contamination
gamma_contamination <- 0.1
rho_contamination <- 0.9
gamma_n <- floor(n*gamma_contamination)

# Scenarios
sample_dim <- rbind(c(50,30), c(100,200), c(100,1000))
colnames(sample_dim) <- c("n", "p")
################################################################################
# Checks

stopifnot(p >= 16)
stopifnot(n%%2 == 0)
################################################################################
# Preparing the simulation data

theta1 <- rep(1,p)
theta2 <- c(rep(3,5), rep(1, p - 5))

tmp_i <- sapply(1:p, function(x) rep(x, p))
tmp <- abs(tmp_i - t(tmp_i))

Sigma <- rho^tmp

rm(tmp, tmp_i)

# Contaminated variables
theta_contaminated <- rep(0.5, p)
theta_contaminated[1:5] <- 3
# Contaminated correlation matrix
Sigma_contaminated <- matrix(rho_contamination, p,p)
diag(Sigma_contaminated) <- 1
################################################################################
beta_tmp <- rep(0,p)
beta_tmp[c(1,3,5,11,13)] <- -0.5
beta_tmp[2] <- 1
beta_tmp[16] = 1.5
beta <- c(-1, beta_tmp)
names(beta) <- coeff_names
################################################################################
# Sample the log normal distribution
set.seed(123)
W1 <- MASS::mvrnorm(n = n/2, mu = theta1, Sigma = Sigma) %>%
  exp()
set.seed(321)
W2 <- MASS::mvrnorm(n = n/2, mu = theta2, Sigma = Sigma) %>%
  exp()
set.seed(234)
W_contaminated <- MASS::mvrnorm(n = gamma_n, mu = theta_contaminated, Sigma = Sigma_contaminated) %>%
  exp()

colnames(W1) <- coeff_names[-1]
colnames(W2) <- coeff_names[-1]
colnames(W_contaminated) <- coeff_names[-1]

W1 <- as_tibble(W1) %>%
  mutate(data_type = "W1")
W2 <- as_tibble(W2) %>%
  mutate(data_type = "W2")
W_contaminated <- as_tibble(W_contaminated) %>%
  mutate(data_type = "W_contaminated")

W <- bind_rows(as_tibble(W1), as_tibble(W2))

# W <- rbind(W1, W2)
# Z <- log(W/(rowSums(W)[1:4]))
# rowSums(W)[1:3]
################################################################################
Z <- W %>%
  # Log transform the data
  mutate(row_sum = W %>%
           dplyr::select(-data_type) %>%
           rowSums(),
         across(starts_with("p_"), ~ log(.x/row_sum))) %>%
  # Remove row sums as not needed anymore
  dplyr::select(-row_sum) %>%
  # Add binary response and intercept
  mutate(y = ifelse(data_type == "W1", 0, 1),
         p_0 = 1,
         id = paste0("clean_",row_number())) %>%
  relocate(p_0)

# z_beta <- as.numeric(as.matrix(Z[,coeff_names]) %*% beta)

Z <- Z %>%
  # Add mu value
  mutate(mu_beta_z = compute_mu_beta_z(beta = beta, z = as.matrix(Z[,coeff_names])))

# Arrange according to mu_beta_z for group W1
Z <- Z %>%
  filter(data_type == "W1") %>%
  arrange(mu_beta_z) %>%
  bind_rows(Z %>%
              filter(data_type == "W2"))

################################################################################
# Contamination schemes
Z_tmp <- W_contaminated %>%
  # Log transform the data
  mutate(row_sum = W_contaminated %>%
           dplyr::select(-data_type) %>%
           rowSums(),
         across(starts_with("p_"), ~ log(.x/row_sum))) %>%
  # Remove row sums as not needed anymore
  dplyr::select(-row_sum) %>%
  # Add binary response and intercept
  mutate(p_0 = 1,
         id = paste0("cont_", row_number())) %>%
  relocate(p_0)

################################################################################

Z_final <- Z_tmp %>%
  bind_cols(Z %>%
              slice(1:gamma_n) %>%
              dplyr::select(y, mu_beta_z)) %>%
  bind_rows(Z %>%
              slice((gamma_n + 1):nrow(Z)))



