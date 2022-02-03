library(tidyverse)
library(MASS)
library(boot)
library(MCMCpack)

generate_parameters <- function(p, rho = 0.2, rho_contamination = 0.9){
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
  
  # Contaminated variables
  theta_contaminated <- rep(0.5, p)
  theta_contaminated[1:5] <- 3
  
  # Contaminated correlation matrix
  Sigma_contaminated <- matrix(rho_contamination, p,p)
  diag(Sigma_contaminated) <- 1
  
  return(list(theta1 = theta1,
              theta2 = theta2,
              Sigma = Sigma,
              theta_contaminated = theta_contaminated,
              Sigma_contaminated = Sigma_contaminated))
}

generate_beta <- function(p){
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

generate_W <- function(n, p, gamma_contamination = 0.1, 
                       theta1, theta2, theta_contaminated, 
                       Sigma, Sigma_contaminated,
                       seed_select = 1){
  # Check parameter dimensions
  if(length(theta1) != p) {stop("1 Parameter vectors do not have the same length")}
  if(length(theta2) != p) {stop("2 Parameter vectors do not have the same length")}
  if(length(theta_contaminated) != p) {stop("3 Parameter vectors do not have the same length")}
  if(dim(Sigma)[1] != p) {stop("4 Parameter vectors do not have the same length")}
  if(dim(Sigma)[2] != p) {stop("5 Parameter vectors do not have the same length")}
  if(dim(Sigma_contaminated)[1] != p) {stop("6  vectors do not have the same length")}
  if(dim(Sigma_contaminated)[2] != p) {stop("7 Parameter vectors do not have the same length")}
  # Check number of samples
  stopifnot(n%%2 == 0)
  
  # Number of contaminated samples
  gamma_n <- floor(n*gamma_contamination)
  
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
  set.seed(n+p + 456 + seed_select)
  W_contaminated <- MASS::mvrnorm(n = gamma_n, mu = theta_contaminated, Sigma = Sigma_contaminated) %>%
    exp()
  
  colnames(W1) <- purrr::map_chr(1:p, ~paste0("p_",.x))
  colnames(W2) <- purrr::map_chr(1:p, ~paste0("p_",.x))
  colnames(W_test1) <- purrr::map_chr(1:p, ~paste0("p_",.x))
  colnames(W_test2) <- purrr::map_chr(1:p, ~paste0("p_",.x))
  colnames(W_contaminated) <- purrr::map_chr(1:p, ~paste0("p_",.x))
  
  W1 <- as_tibble(W1) %>%
    mutate(data_type = "W1")
  W2 <- as_tibble(W2) %>%
    mutate(data_type = "W2")
  W_test1 <- as_tibble(W_test1) %>%
    mutate(data_type = "W_test1")
  W_test2 <- as_tibble(W_test2) %>%
    mutate(data_type = "W_test2")
  W_contaminated <- as_tibble(W_contaminated) %>%
    mutate(data_type = "W_contaminated")
  
  W <- bind_rows(as_tibble(W1), as_tibble(W2))
  W_test <- bind_rows(as_tibble(W_test1), as_tibble(W_test2))
  
  return(list(W = W, W_test = W_test, W_contaminated = W_contaminated))
}

generate_Z <- function(p, beta, W, W_test, W_contaminated){
  if(length(beta) != (p+1)) {stop("1 beta vector does not have the same length as p")}
  
  coeff_names <- purrr::map_chr(0:p, ~paste0("p_",.x))
  gamma_n <- nrow(W_contaminated)
  
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
    arrange(mu_beta_z) %>%
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
    arrange(mu_beta_z) %>%
    bind_rows(Z %>%
                filter(data_type == "W2")) %>%
    relocate(data_type, id, y, row_sum, mu_beta_z, p_0)
  
  # Contamination schemes
  Z_tmp <- W_contaminated %>%
    # Log transform the data
    # mutate(row_sum = W_contaminated %>%
    #          dplyr::select(-data_type) %>%
    #          rowSums(),
    #        across(starts_with("p_"), ~ log(.x/row_sum))) %>%
    mutate(row_sum = Z %>% 
             slice(1:gamma_n) %>% 
             pull(row_sum),
           across(starts_with("p_"), ~ log(.x/row_sum))) %>%
    # dplyr::select(-row_sum) %>%
    # Add binary response and intercept
    mutate(p_0 = 1,
           id = paste0("cont_", row_number()))
  
  # Replace first gamma_n covariates with contaminated covariates - Then append the rest
  Z_contaminated <- Z_tmp %>%
    bind_cols(Z %>%
                slice(1:gamma_n) %>%
                dplyr::select(y, mu_beta_z)) %>%
    bind_rows(Z %>%
                slice((gamma_n + 1):nrow(Z))) %>%
    relocate(data_type, id, y, row_sum, mu_beta_z, p_0)

  return(list(Z = Z, Z_test = Z_test, Z_contaminated = Z_contaminated))
}

################################################################################
################################################################################
################################################################################

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

################################################################################
################################################################################
################################################################################


generate_parameters3 <- function(n, p, gamma_n, seed){
  stopifnot(p >= 30)
  
  par_5 <- 10
  
  alpha0 <- c(rep(60, 3), rep(120, 4), rep(5, par_5), rep(0.1, p - (3 + 4 + par_5)))
  alpha1 <- c(rep(0.1, 2), rep(90, 5), rep(60, 5), rep(5, par_5), rep(0.1, p - (5 + 5 + par_5 + 2)))
  alpha_contaminated <- c(rep(70, 4), rep(85,4), rep(5, par_5), rep(0.1, p - (4+4 + par_5)))
  
  set.seed(seed + 1323)
  p_dirichlet0 <- MCMCpack::rdirichlet(n/2, alpha0)
  set.seed(seed + 2334)
  p_dirichlet1 <- MCMCpack::rdirichlet(n/2, alpha1)
  set.seed(seed + 1345)
  p_dirichlet_contaminated <- MCMCpack::rdirichlet(gamma_n, alpha_contaminated)
  
  return(list(
    p_dirichlet0 = p_dirichlet0,
    p_dirichlet1 = p_dirichlet1,
    p_dirichlet_contaminated = p_dirichlet_contaminated
  ))
}

generate_multinomial <- function(n, size, prob, seed =1){
  set.seed(seed)
  X <- t(rmultinom(n = n, size = size, prob = prob))
  colnames(X) = purrr::map_chr(1:ncol(X), ~paste0("p_",.x))
  
  return(as_tibble(X))
}

generate_W3 <- function(n, gamma_n, p_dirichlet0, p_dirichlet1, p_dirichlet_contaminated){
  mu <- 6000
  sigma <- 600
  
  set.seed(4444)
  n_samples_train0 <- round(rnorm(n/2, mu, sigma))
  set.seed(3452)
  n_samples_train1 <- round(rnorm(n/2, mu, sigma))
  set.seed(4443)
  n_samples_test0 <- round(rnorm(n/2, mu, sigma))
  set.seed(3453)
  n_samples_test1 <- round(rnorm(n/2, mu, sigma))
  set.seed(4442)
  n_samples_contaminated <- round(rnorm(gamma_n, sigma))
  
  W0 <- purrr::map_df(1:(n/2), ~generate_multinomial(n = 1, size = n_samples_train0[.x], prob = p_dirichlet0[.x,], seed = n_samples_train0[.x])) %>%
    mutate(across(everything(), ~ifelse(.x == 0, 0.5, .x))) %>%
    normalize_rows() %>%
    mutate(data_type = "W0")
  W1 <- purrr::map_df(1:(n/2), ~generate_multinomial(n = 1, size = n_samples_train1[.x], prob = p_dirichlet1[.x,], seed = n_samples_train1[.x])) %>%
    mutate(across(everything(), ~ifelse(.x == 0, 0.5, .x))) %>%
    normalize_rows() %>%
    mutate(data_type = "W1")
  
  W_test0 <- purrr::map_df(1:(n/2), ~generate_multinomial(n = 1, size = n_samples_test0[.x], prob = p_dirichlet0[.x,], seed = n_samples_test0[.x])) %>%
    mutate(across(everything(), ~ifelse(.x == 0, 0.5, .x))) %>%
    normalize_rows() %>%
    mutate(data_type = "W_test0")
  W_test1 <- purrr::map_df(1:(n/2), ~generate_multinomial(n = 1, size = n_samples_test1[.x], prob = p_dirichlet1[.x,], seed = n_samples_test1[.x])) %>%
    mutate(across(everything(), ~ifelse(.x == 0, 0.5, .x))) %>%
    normalize_rows() %>%
    mutate(data_type = "W_test1")
  
  W_contaminated <- purrr::map_df(1:gamma_n, ~generate_multinomial(n = 1, size = n_samples_contaminated[.x], prob = p_dirichlet_contaminated[.x,], seed = n_samples_contaminated[.x])) %>%
    mutate(across(everything(), ~ifelse(.x == 0, 0.5, .x))) %>%
    normalize_rows() %>%
    mutate(data_type = "W_contaminated")
  
  W <- bind_rows(W0, W1)
  W_test <- bind_rows(W_test0, W_test1)
  
  return(list(W = W, W_test = W_test, W_contaminated = W_contaminated))
}



################################################################################
################################################################################
generate_Z3 <- function(p, gamma_n, beta, W, W_test, W_contaminated){
  if(length(beta) != (p+1)) {stop("1 beta vector does not have the same length as p")}
  
  coeff_names <- purrr::map_chr(0:p, ~paste0("p_",.x))
  
  Z_test <- W_test %>%
    mutate(row_sum = W_test %>%
             dplyr::select(-data_type) %>%
             rowSums(),
           across(starts_with("p_"), ~ log(.x/row_sum))) %>%
    mutate(y = ifelse(data_type == "W_test0", 0, 1),
           p_0 = 1,
           id = paste0("clean_", row_number())) %>%
    relocate(p_0)
  
  # Add mu value
  Z_test <- Z_test %>%
    mutate(mu_beta_z = compute_mu_beta_z(beta = beta, z = as.matrix(Z_test[,coeff_names])))
  
  # Arrange according to mu_beta_z for group W1
  Z_test <- Z_test %>%
    filter(data_type == "W_test0") %>%
    arrange(mu_beta_z) %>%
    bind_rows(Z_test %>%
                filter(data_type == "W_test1")) %>%
    relocate(data_type, id, y, row_sum, mu_beta_z, p_0)
  
  Z <- W %>%
    mutate(row_sum = W_test %>%
             dplyr::select(-data_type) %>%
             rowSums(),
           across(starts_with("p_"), ~ log(.x/row_sum))) %>%
    mutate(y = ifelse(data_type == "W0", 0, 1),
           p_0 = 1,
           id = paste0("clean_", row_number())) %>%
    relocate(p_0)
  
  # Add mu value
  Z <- Z %>%
    mutate(mu_beta_z = compute_mu_beta_z(beta = beta, z = as.matrix(Z[,coeff_names])))
  
  # Arrange according to mu_beta_z for group W1
  Z <- Z %>%
    filter(data_type == "W0") %>%
    arrange(mu_beta_z) %>%
    bind_rows(Z %>%
                filter(data_type == "W1")) %>%
    relocate(data_type, id, y, row_sum, mu_beta_z, p_0)
  
  # Contamination schemes
  Z_tmp <- W_contaminated %>%
    # Log transform the data
    # mutate(row_sum = W_contaminated %>%
    #          dplyr::select(-data_type) %>%
    #          rowSums(),
    #        across(starts_with("p_"), ~ log(.x/row_sum))) %>%
    mutate(row_sum = Z %>% 
             slice(1:gamma_n) %>% 
             pull(row_sum),
           across(starts_with("p_"), ~ log(.x/row_sum))) %>%
    # dplyr::select(-row_sum) %>%
    # Add binary response and intercept
    mutate(p_0 = 1,
           id = paste0("cont_", row_number()))
  
  # Replace first gamma_n covariates with contaminated covariates - Then append the rest
  Z_contaminated <- Z_tmp %>%
    bind_cols(Z %>%
                slice(1:gamma_n) %>%
                dplyr::select(y, mu_beta_z)) %>%
    bind_rows(Z %>%
                slice((gamma_n + 1):nrow(Z))) %>%
    relocate(data_type, id, y, row_sum, mu_beta_z, p_0)
  
  return(list(Z = Z, Z_test = Z_test, Z_contaminated = Z_contaminated))
}

