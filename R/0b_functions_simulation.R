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
    mutate(row_sum = W_contaminated %>%
             dplyr::select(-data_type) %>%
             rowSums(),
           across(starts_with("p_"), ~ log(.x/row_sum))) %>%
    # mutate(row_sum = Z %>% 
    #          slice(1:gamma_n) %>% 
    #          pull(row_sum),
    #        across(starts_with("p_"), ~ log(.x/row_sum))) %>%
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
# Functions for own simulation study
################################################################################
mix_geom_pois <- function(n, par_geom, par_pois, p1 = 0.9, seed){
  set.seed(seed + 1000)
  alpha <- as.numeric(runif(n = n) < p1)
  
  set.seed(seed + 2000)
  return(alpha*rgeom(n = n, prob = 1/par_geom) + (1-alpha)*rpois(n = n, lambda = par_pois))
}

generate_W_microbiome <- function(n, p, gamma_n, par_0, par_contaminated, seed){
  coeff_names <- paste0("p_", seq(1:p))
  
  par_train_0 <- par_0 %>%
    mutate(n = n,
           seed = seed*1000+150 + row_number())
  par_test_0 <- par_0 %>%
    mutate(n = n,
           seed = seed*1000+450 + row_number())
  par_contaminated <- par_contaminated %>%
    mutate(n = gamma_n,
           seed = seed*1000+750 + row_number())
  
  W0 <- suppressMessages(purrr::pmap_dfc(par_train_0 %>% dplyr::select(n, par_geom, par_pois, seed), mix_geom_pois))
  W_test0 <- suppressMessages(purrr::pmap_dfc(par_test_0 %>% dplyr::select(n, par_geom, par_pois, seed), mix_geom_pois))
  W_contaminated <- suppressMessages(purrr::pmap_dfc(par_contaminated %>% dplyr::select(n, par_geom, par_pois, seed), mix_geom_pois))
  
  colnames(W0) <- coeff_names
  colnames(W_test0) <- coeff_names
  colnames(W_contaminated) <- coeff_names
  
  W <- W0 %>%
    mutate(across(everything(), ~ifelse(.x == 0, 0.5, .x)),
           data_type = "W0") %>%
    relocate(data_type)
  W_test <- W_test0 %>%
    mutate(across(everything(), ~ifelse(.x == 0, 0.5, .x)),
           data_type = "W_test0") %>%
    relocate(data_type)
  W_contaminated <- W_contaminated %>%
    mutate(across(everything(), ~ifelse(.x == 0, 0.5, .x)),
           data_type = "W_contaminated") %>%
    relocate(data_type)

  # par_train_1 <- par_1 %>%
  #   mutate(n = n,
  #          seed = seed*1000+300 + row_number())
  # par_test_1 <- par_1 %>%
  #   mutate(n = n,
  #          seed = seed*1000+600 + row_number())
  
  # W1 <- purrr::pmap_dfc(par_train_1 %>% dplyr::select(n, par_geom, par_pois, seed), mix_geom_pois)
  # W_test1 <- purrr::pmap_dfc(par_test_1 %>% dplyr::select(n, par_geom, par_pois, seed), mix_geom_pois)
  # colnames(W1) <- coeff_names
  # colnames(W_test1) <- coeff_names
  
  # W1 <- W1 %>%
  #   mutate(across(everything(), ~ifelse(.x == 0, 0.5, .x))) %>%
  #   # normalize_rows() %>%
  #   mutate(data_type = "W1")
  # W_test1 <- W_test1 %>%
  #   mutate(across(everything(), ~ifelse(.x == 0, 0.5, .x))) %>%
  #   # normalize_rows() %>%
  #   mutate(data_type = "W_test1")
  
  return(list(W = W, W_test = W_test, W_contaminated = W_contaminated))
}

generate_Z_microbiome <- function(n, p, gamma_n, beta, W, W_test, W_contaminated, seed){
  if(length(beta) != (p+1)) {stop("1 beta vector does not have the same length as p")}
  
  coeff_names <- purrr::map_chr(0:p, ~paste0("p_",.x))
  
  #########################################################
  
  Z_test <- W_test %>%
    mutate(row_sum = W_test %>%
             dplyr::select(starts_with("p_")) %>%
             rowSums(),
           across(starts_with("p_"), ~ log(.x/row_sum))) %>%
    mutate(p_0 = 1,
           id = paste0("clean_", row_number()))
  
  set.seed(n+p + 123 + seed)
  # Add mu value
  Z_test <- Z_test %>%
    mutate(mu_beta_z = compute_mu_beta_z(beta = beta, z = as.matrix(Z_test[,coeff_names])),
           y = as.numeric(runif(nrow(Z_test)) < mu_beta_z)) %>%
    arrange(mu_beta_z) %>%
    relocate(data_type, id, y, row_sum, mu_beta_z, p_0)
  
  #########################################################
  
  Z <- W %>%
    mutate(row_sum = W %>%
             dplyr::select(starts_with("p_")) %>%
             rowSums(),
           across(starts_with("p_"), ~ log(.x/row_sum))) %>%
    mutate(p_0 = 1,
           id = paste0("clean_", row_number()))
  
  set.seed(n+p + 567 + seed)
  # Add mu value
  Z <- Z %>%
    mutate(mu_beta_z = compute_mu_beta_z(beta = beta, z = as.matrix(Z[,coeff_names])),
           y = as.numeric(runif(nrow(Z)) < mu_beta_z)) %>%
    arrange(mu_beta_z) %>%
    relocate(data_type, id, y, row_sum, mu_beta_z, p_0)
  
  #########################################################
  
  # Contamination schemes
  Z_tmp <- W_contaminated %>%
    # Log transform the data
    mutate(row_sum = W_contaminated %>%
             dplyr::select(starts_with("p_")) %>%
             rowSums(),
           across(starts_with("p_"), ~ log(.x/row_sum))) %>%
    mutate(p_0 = 1,
           id = paste0("cont_", row_number()))
  
  # Add mu value
  Z_tmp <- Z_tmp %>%
    mutate(mu_beta_z = compute_mu_beta_z(beta = beta, z = as.matrix(Z_tmp[,coeff_names])))
  
  # Replace first gamma_n covariates with contaminated covariates - Then append the rest
  Z_contaminated <- Z_tmp %>%
    bind_cols(Z %>%
                slice(1:gamma_n) %>%
                dplyr::select(y)) %>%
    bind_rows(Z %>%
                slice((gamma_n + 1):nrow(Z))) %>%
    relocate(data_type, id, y, row_sum, mu_beta_z, p_0)
  
  return(list(Z = Z, Z_test = Z_test, Z_contaminated = Z_contaminated))
}

