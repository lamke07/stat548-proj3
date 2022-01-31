################################################################################
rm(list = ls())

source("0a_functions_data.R")
source("0b_functions_utils.R")

dir.create("data")
################################################################################
#' Example Hyperparameters
# n <- 20
# p <- 17
# rho <- 0.2
# rho_contamination <- 0.9
# gamma_contamination <- 0.1

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
  
  sim_par <- generate_parameters(p)
  sim_beta <- generate_beta(p)
  sim_W <- generate_W(n = n, p = p,
                      theta1 = sim_par$theta1, theta2 = sim_par$theta2, theta_contaminated = sim_par$theta_contaminated,
                      Sigma = sim_par$Sigma, Sigma_contaminated = sim_par$Sigma_contaminated,
                      seed_select = i)
  sim_Z <- generate_Z(p = p, beta = sim_beta$beta, W = sim_W$W, W_test = sim_W$W_test, W_contaminated = sim_W$W_contaminated)
  
  sim_all <- c(sim_par, sim_beta, sim_W, sim_Z)
  return(sim_all)
})

sim2 <- lapply(1:100, function(i){
  n = sample_dim[2, "n"]
  p = sample_dim[2, "p"]
  
  cat("\nSimulating for parameters (n,p) = (", n, ", ", p, ")\n")
  
  sim_par <- generate_parameters(p)
  sim_beta <- generate_beta(p)
  sim_W <- generate_W(n = n, p = p,
                      theta1 = sim_par$theta1, theta2 = sim_par$theta2, theta_contaminated = sim_par$theta_contaminated,
                      Sigma = sim_par$Sigma, Sigma_contaminated = sim_par$Sigma_contaminated,
                      seed_select = i*10)
  sim_Z <- generate_Z(p = p, beta = sim_beta$beta, W = sim_W$W, W_test = sim_W$W_test, W_contaminated = sim_W$W_contaminated)
  
  sim_all <- c(sim_par, sim_beta, sim_W, sim_Z)
  return(sim_all)
})

# sim3 <- lapply(1:100, function(i){
#   n = sample_dim[3, "n"]
#   p = sample_dim[3, "p"]
#   
#   cat("\nSimulating for parameters (n,p) = (", n, ", ", p, ")\n")
#   
#   sim_par <- generate_parameters(p)
#   sim_beta <- generate_beta(p)
#   sim_W <- generate_W(n = n, p = p,
#                       theta1 = sim_par$theta1, theta2 = sim_par$theta2, theta_contaminated = sim_par$theta_contaminated,
#                       Sigma = sim_par$Sigma, Sigma_contaminated = sim_par$Sigma_contaminated,
#                       seed_select = i*100)
#   sim_Z <- generate_Z(p = p, beta = sim_beta$beta, W = sim_W$W, W_contaminated = sim_W$W_contaminated)
#   
#   sim_all <- c(sim_par, sim_beta, sim_W, sim_Z)
#   return(sim_all)
# })

saveRDS(sim1, "data/sim1.RDS")
saveRDS(sim2, "data/sim2.RDS")
# saveRDS(sim3, "data/sim3.RDS")
