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

for(i in (1:nrow(sample_dim))){
  n = sample_dim[i, "n"]
  p = sample_dim[i, "p"]
  
  cat("\nSimulating for parameters (n,p) = (", n, ", ", p, ")\n")
  
  sim_par <- generate_parameters(p)
  sim_beta <- generate_beta(p)
  sim_W <- generate_W(n = n, p = p,
                      theta1 = sim_par$theta1, theta2 = sim_par$theta2, theta_contaminated = sim_par$theta_contaminated,
                      Sigma = sim_par$Sigma, Sigma_contaminated = sim_par$Sigma_contaminated)
  sim_Z <- generate_Z(p = p, beta = sim_beta$beta, W = sim_W$W, W_contaminated = sim_W$W_contaminated)
  
  sim_all <- c(sim_par, sim_beta, sim_W, sim_Z)
  saveRDS(sim_all, file = paste0("data/", i, "_sim_all.RDS"))
  readr::write_csv(sim_Z$Z, paste0("data/", i, "_sim_Z.csv"))
  readr::write_csv(sim_Z$Z_contaminated, paste0("data/", i, "_sim_Z_contaminated.csv"))
}