library(tidyverse)
library(MASS)
library(boot)
library(MCMCpack)

generate_parameters2 <- function(p, rho = 0.2){
  # Checks for dimension
  stopifnot(p >= 16)
  
  # means of the lognormal distribution
  theta1 <- rep(1,p)
  theta2 <- c(rep(3,5), rep(1, p - 5))
  
  # Covariance matrix of the lognormal distribution
  tmp_i <- sapply(1:p, function(x) rep(x, p))
  tmp <- abs(tmp_i - t(tmp_i))
  Sigma <- rho^tmp
  rm(tmp, tmp_i)
  
  return(list(theta1 = theta1,
              theta2 = theta2,
              Sigma = Sigma))
}

generate_beta2 <- function(p){
  # Checks for dimension
  stopifnot(p >= 16)
  
  beta_tmp <- rep(0,p)
  beta_tmp[c(1,3,5,11,13)] <- -0.5
  beta_tmp[2] <- 1
  beta_tmp[16] <- 1.5
  beta <- c(-1, beta_tmp)
  names(beta) <- purrr::map_chr(0:p, ~paste0("p_",.x))
  
  return(list(beta = beta))
}

generate_W2 <- function(n, p,
                        theta1, theta2, Sigma,
                        seed_select = 1){
  # Check parameter dimensions
  if(length(theta1) != p) {stop("1 Parameter vectors do not have the same length")}
  if(length(theta2) != p) {stop("2 Parameter vectors do not have the same length")}
  if(dim(Sigma)[1] != p) {stop("4 Parameter vectors do not have the same length")}
  if(dim(Sigma)[2] != p) {stop("5 Parameter vectors do not have the same length")}
  # Check number of samples
  stopifnot(n%%2 == 0)
  
  # Sample the log normal distribution
  set.seed(n+p + 123 + seed_select)
  W1 <- MASS::mvrnorm(n = n/2, mu = theta1, Sigma = Sigma) %>%
    exp()
  set.seed(n+p + 231 + seed_select)
  W2 <- MASS::mvrnorm(n = n/2, mu = theta2, Sigma = Sigma) %>%
    exp()
  set.seed(n+p + 567 + seed_select)
  W_test1 <- MASS::mvrnorm(n = n/2, mu = theta1, Sigma = Sigma) %>%
    exp()
  set.seed(n+p + 987 + seed_select)
  W_test2 <- MASS::mvrnorm(n = n/2, mu = theta2, Sigma = Sigma) %>%
    exp()
  
  colnames(W1) <- purrr::map_chr(1:p, ~paste0("p_",.x))
  colnames(W2) <- purrr::map_chr(1:p, ~paste0("p_",.x))
  colnames(W_test1) <- purrr::map_chr(1:p, ~paste0("p_",.x))
  colnames(W_test2) <- purrr::map_chr(1:p, ~paste0("p_",.x))
  
  W1 <- as_tibble(W1) %>%
    mutate(data_type = "W1")
  W2 <- as_tibble(W2) %>%
    mutate(data_type = "W2")
  W_test1 <- as_tibble(W_test1) %>%
    mutate(data_type = "W_test1")
  W_test2 <- as_tibble(W_test2) %>%
    mutate(data_type = "W_test2")
  
  W <- bind_rows(as_tibble(W1), as_tibble(W2))
  W_test <- bind_rows(as_tibble(W_test1), as_tibble(W_test2))
  
  return(list(W = W, W_test = W_test))
}


generate_X_dirichlet <- function(n, p, seed_select, gamma_contamination = 0.2){
  stopifnot(p >= 11)
  
  gamma_n <- floor(n*gamma_contamination)
  set.seed(5555)
  alpha <- rep(1, p) + c(rep(0, 10), sample(1:10, size = p-10, replace = TRUE))
  
  set.seed(n + p + 543 + seed_select)
  X_dirichlet <- MCMCpack::rdirichlet(gamma_n, alpha)
  colnames(X_dirichlet) <- purrr::map_chr(1:p, ~paste0("p_",.x))
  
  X_dirichlet <- X_dirichlet %>%
    as_tibble() %>%
    mutate(data_type = "W_dirichlet")
  
  return(list(X_dirichlet = X_dirichlet))
}


generate_Z2 <- function(p, beta, W, W_test, X_dirichlet){
  if(length(beta) != (p+1)) {stop("1 beta vector does not have the same length as p")}
  
  gamma_n <- nrow(X_dirichlet)
  coeff_names <- purrr::map_chr(0:p, ~paste0("p_",.x))
  
  # Log transform the test data
  Z_test <- W_test %>%
    mutate(row_sum = W_test %>%
             dplyr::select(-data_type) %>%
             rowSums(),
           across(starts_with("p_"), ~ log(.x/row_sum))) %>%
    # dplyr::select(-row_sum) %>%
    # Add binary response and intercept
    mutate(y = ifelse(data_type == "W_test1", 0, 1),
           p_0 = 1,
           id = paste0("clean_",row_number())) %>%
    relocate(p_0)
  
  # Add mu value
  Z_test <- Z_test %>%
    mutate(mu_beta_z = compute_mu_beta_z(beta = beta, z = as.matrix(Z_test[,coeff_names])))
  
  # Arrange according to mu_beta_z for group W1
  Z_test <- Z_test %>%
    filter(data_type == "W_test1") %>%
    arrange(-mu_beta_z) %>%
    bind_rows(Z_test %>%
                filter(data_type == "W_test2")) %>%
    relocate(data_type, id, y, row_sum, mu_beta_z, p_0)
  
  # Log transform the training data
  Z <- W %>%
    mutate(row_sum = W %>%
             dplyr::select(-data_type) %>%
             rowSums(),
           across(starts_with("p_"), ~ log(.x/row_sum))) %>%
    # dplyr::select(-row_sum) %>%
    # Add binary response and intercept
    mutate(y = ifelse(data_type == "W1", 0, 1),
           p_0 = 1,
           id = paste0("clean_",row_number()))
  
  # Add mu value
  Z <- Z %>%
    mutate(mu_beta_z = compute_mu_beta_z(beta = beta, z = as.matrix(Z[,coeff_names])))
  
  # Arrange according to mu_beta_z for group W1
  Z <- Z %>%
    filter(data_type == "W1") %>%
    arrange(-mu_beta_z) %>%
    bind_rows(Z %>%
                filter(data_type == "W2")) %>%
    relocate(data_type, id, y, row_sum, mu_beta_z, p_0)
  
  # Replace first gamma_n covariates with contaminated covariates - Then append the rest
  Z_contaminated <- X_dirichlet %>%
    mutate(p_0 = 1,
           id = paste0("dir_", row_number())) %>%
    bind_cols(Z %>%
                slice(1:gamma_n) %>%
                dplyr::select(y, mu_beta_z)) %>%
    bind_rows(Z %>%
                slice((gamma_n + 1):nrow(Z))) %>%
    relocate(data_type, id, y, row_sum, mu_beta_z, p_0)
  
  return(list(Z = Z, Z_test = Z_test, Z_contaminated = Z_contaminated))
}

