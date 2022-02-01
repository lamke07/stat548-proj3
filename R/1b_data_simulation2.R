################################################################################
rm(list = ls())

source("0a_functions_utils.R")
source("0c_functions_simulation2.R")

dir.create("data")
################################################################################
#' Example Hyperparameters
# n <- 20
# p <- 17
# rho <- 0.2
# gamma_contamination <- 0.2

################################################################################
# Scenarios
sample_dim <- rbind(c(50,30),
                    c(100,200),
                    c(100,1000))
colnames(sample_dim) <- c("n", "p")
# i = 1
sim1 <- lapply(1:100, function(i){
  n = sample_dim[1, "n"]
  p = sample_dim[1, "p"]
  
  cat("\nSimulating for parameters (n,p) = (", n, ", ", p, ")\n")
  
  sim_par <- generate_parameters2(p)
  sim_beta <- generate_beta2(p)
  sim_W <- generate_W2(n = n, p = p,
                       theta1 = sim_par$theta1, theta2 = sim_par$theta2,
                       Sigma = sim_par$Sigma,
                       seed_select = i+100)
  sim_X_dirichlet <- generate_X_dirichlet(n = n, p = p, seed_select = i + 6)
  sim_Z <- generate_Z2(p = p, beta = sim_beta$beta, W = sim_W$W, W_test = sim_W$W_test, X_dirichlet = sim_X_dirichlet$X_dirichlet)
  
  sim_all <- c(sim_par, sim_beta, sim_W, sim_X_dirichlet, sim_Z)
  return(sim_all)
})

sim2 <- lapply(1:100, function(i){
  n = sample_dim[2, "n"]
  p = sample_dim[2, "p"]
  
  cat("\nSimulating for parameters (n,p) = (", n, ", ", p, ")\n")
  
  sim_par <- generate_parameters2(p)
  sim_beta <- generate_beta2(p)
  sim_W <- generate_W2(n = n, p = p,
                       theta1 = sim_par$theta1, theta2 = sim_par$theta2,
                       Sigma = sim_par$Sigma,
                       seed_select = i+100)
  sim_X_dirichlet <- generate_X_dirichlet(n = n, p = p, seed_select = i + 4)
  sim_Z <- generate_Z2(p = p, beta = sim_beta$beta, W = sim_W$W, W_test = sim_W$W_test, X_dirichlet = sim_X_dirichlet$X_dirichlet)
  
  sim_all <- c(sim_par, sim_beta, sim_W, sim_X_dirichlet, sim_Z)
  return(sim_all)
})

# sim3 <- lapply(1:100, function(i){
#   n = sample_dim[3, "n"]
#   p = sample_dim[3, "p"]
#   
#   cat("\nSimulating for parameters (n,p) = (", n, ", ", p, ")\n")
#   
#   sim_par <- generate_parameters2(p)
#   sim_beta <- generate_beta2(p)
#   sim_W <- generate_W2(n = n, p = p,
#                        theta1 = sim_par$theta1, theta2 = sim_par$theta2,
#                        Sigma = sim_par$Sigma,
#                        seed_select = i+100)
#   sim_X_dirichlet <- generate_X_dirichlet(n = n, p = p, seed_select = i + 5)
#   sim_Z <- generate_Z2(p = p, beta = sim_beta$beta, W = sim_W$W, W_test = sim_W$W_test, X_dirichlet = sim_X_dirichlet$X_dirichlet)
#   
#   sim_all <- c(sim_par, sim_beta, sim_W, sim_X_dirichlet, sim_Z)
#   return(sim_all)
# })

saveRDS(sim1, "data/extended_sim1.RDS")
saveRDS(sim2, "data/extended_sim2.RDS")
# saveRDS(sim3, "data/extended_sim3.RDS")
