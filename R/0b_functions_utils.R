library(tidyverse)
library(MASS)
library(boot)

compute_mu_beta_z <- function(beta, z){
  # Input: vector[p] beta, matrix[n,p] z
  # Check all dimensions are the same
  stopifnot(dim(z)[2] == length(beta))
  boot::inv.logit(as.numeric(z %*% beta))
}