################################################################################
rm(list = ls())
library(MCMCpack)

source("0a_functions_utils.R")
source("0b_functions_simulation.R")

################################################################################
# Scenarios
sample_dim <- rbind(c(50,30),
                    c(100,150))
colnames(sample_dim) <- c("n", "p")
# i = 1
sim1 <- lapply(1:100, function(i){
  n = sample_dim[1, "n"]
  p = sample_dim[1, "p"]
  gamma_n <- floor(n*0.1)
  
  cat("\nSimulating for parameters (n,p) = (", n, ", ", p, ")\n")
  
  sim_par <- generate_parameters3(n = n, p = p, gamma_n = gamma_n, seed = i)
  sim_beta <- generate_beta(p)
  sim_W <- generate_W3(n = n, gamma_n = gamma_n,
                       p_dirichlet0 = sim_par$p_dirichlet0, p_dirichlet1 = sim_par$p_dirichlet1,
                       p_dirichlet_contaminated = sim_par$p_dirichlet_contaminated)
  sim_Z <- generate_Z3(p = p, gamma_n = gamma_n, beta = sim_beta$beta, W = sim_W$W, W_test = sim_W$W_test, W_contaminated = sim_W$W_contaminated)
  
  sim_all <- c(sim_par, sim_beta, sim_W, sim_Z)
  return(sim_all)
})


sim2 <- lapply(1:100, function(i){
  n = sample_dim[2, "n"]
  p = sample_dim[2, "p"]
  gamma_n <- floor(n*0.1)
  
  cat("\nSimulating for parameters (n,p) = (", n, ", ", p, ")\n")
  
  sim_par <- generate_parameters3(n = n, p = p, gamma_n = gamma_n, seed = i)
  sim_beta <- generate_beta(p)
  sim_W <- generate_W3(n = n, gamma_n = gamma_n,
                       p_dirichlet0 = sim_par$p_dirichlet0, p_dirichlet1 = sim_par$p_dirichlet1,
                       p_dirichlet_contaminated = sim_par$p_dirichlet_contaminated)
  sim_Z <- generate_Z3(p = p, gamma_n = gamma_n, beta = sim_beta$beta, W = sim_W$W, W_test = sim_W$W_test, W_contaminated = sim_W$W_contaminated)
  
  sim_all <- c(sim_par, sim_beta, sim_W, sim_Z)
  return(sim_all)
})

saveRDS(sim1, "data/extended2_sim1.RDS")
saveRDS(sim2, "data/extended2_sim2.RDS")

################################################################################
################################################################################
################################################################################
################################################################################
################################################################################


################################################################################





