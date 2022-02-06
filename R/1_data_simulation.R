################################################################################
rm(list = ls())

source("0a_functions_utils.R")
source("0b_functions_simulation.R")

dir.create("data")

################################################################################
# Replicate Simulation study
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
  
  cat("\nSimulating for parameters (n,p) = (", n, ", ", p, ")")
  
  sim_par <- generate_parameters(p)
  sim_beta <- generate_beta(p)
  sim_W <- generate_W(n = n, p = p,
                      theta1 = sim_par$theta1, theta2 = sim_par$theta2, theta_contaminated = sim_par$theta_contaminated,
                      Sigma = sim_par$Sigma, Sigma_contaminated = sim_par$Sigma_contaminated,
                      seed_select = 100+i)
  sim_Z <- generate_Z(p = p, beta = sim_beta$beta, W = sim_W$W, W_test = sim_W$W_test, W_contaminated = sim_W$W_contaminated)
  
  sim_all <- c(sim_par, sim_beta, sim_W, sim_Z)
  return(sim_all)
})

sim2 <- lapply(1:100, function(i){
  n = sample_dim[2, "n"]
  p = sample_dim[2, "p"]
  
  cat("\nSimulating for parameters (n,p) = (", n, ", ", p, ")")
  
  sim_par <- generate_parameters(p)
  sim_beta <- generate_beta(p)
  sim_W <- generate_W(n = n, p = p,
                      theta1 = sim_par$theta1, theta2 = sim_par$theta2, theta_contaminated = sim_par$theta_contaminated,
                      Sigma = sim_par$Sigma, Sigma_contaminated = sim_par$Sigma_contaminated,
                      seed_select = 200 + i)
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

################################################################################
# Own Simulation study
print("Starting own simulation study...")
################################################################################
# Create dataframe with cases for geometric and poisson parameters
par_all <- rbind(c(15, 100),
                 c(50, 100),
                 c(50, 300),
                 c(300, 500),
                 c(300, 1000),
                 c(500, 5000))
colnames(par_all) <- c("par_geom", "par_pois")
par_all <- par_all %>% as_tibble() %>%
  mutate(case_id = paste0("par_case_", seq(1,6))) %>%
  relocate(case_id)

plt_df <- suppressMessages(purrr::map2_dfc(par_all$par_geom, par_all$par_pois, ~mix_geom_pois(n = 1000, par_geom = .x, par_pois = .y, seed = 9)))
colnames(plt_df) <- paste0("Case ", seq(1,6))
p <- ggplot(plt_df %>%
              pivot_longer(everything(), names_to = "case", values_to = "samples")) +
  geom_histogram(aes(x = samples)) +
  facet_wrap(~case, scale = "free") +
  labs(x = "Samples", y = "Count", title = "Sampled Distributions") +
  theme_light()
ggsave(p, filename = "figures/sampled_cases.pdf", width = 8, height = 5, dpi = "retina")



# Select distribution settings for case 0 and case 1
case_0 <- c(rep(1, 5), rep(2, 10), rep(3, 25), rep(4, 5), rep(5, 5), rep(6, 10))
case_contaminated <- c(rep(3, 20), rep(2, 10), rep(5, 10), rep(1, 5), rep(6, 15))

par_0 <- paste0("par_case_", case_0) %>%
  as_tibble() %>%
  rename(case_id = value) %>%
  left_join(par_all, by = "case_id")

par_contaminated <- paste0("par_case_", case_contaminated) %>%
  as_tibble() %>%
  rename(case_id = value) %>%
  left_join(par_all, by = "case_id")

################################################################################
n = 150
p <- nrow(par_0)
gamma_n <- floor(n*0.1)

# i = 100
sim_own <- lapply(1:100, function(i){
  
  cat("\nSimulating for parameters (n,p) = (", n, ", ", p, ")")
  
  sim_beta <- generate_beta(p)
  sim_W <- generate_W_microbiome(n = n, p = p, gamma_n = gamma_n, par_0 = par_0, par_contaminated = par_contaminated, seed = i)
  sim_Z <- generate_Z_microbiome(n = n, p = p, gamma_n = gamma_n, beta = sim_beta$beta, W = sim_W$W, W_test = sim_W$W_test, W_contaminated = sim_W$W_contaminated, seed = i*1000)
  
  sim_all <- c(sim_beta, sim_W, sim_Z)
  return(sim_all)
})
saveRDS(sim_own, "data/sim_own.RDS")

################################################################################
# Real data set
print("Saving real data sets")
################################################################################
# Script to save compositional datasets

readr::write_csv(selbal::HIV, "data/selbal_HIV.csv")
# readr::write_csv(selbal::Crohn, "data/selbal_crohn.csv")
# readr::write_csv(as_tibble(selbal::sCD14), "data/selbal_sCD14.csv")
