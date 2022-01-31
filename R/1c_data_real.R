# Script to save compositional datasets
library(selbal)
library(dplyr)
library(readr)

readr::write_csv(selbal::Crohn, "data/selbal_crohn.csv")
readr::write_csv(selbal::HIV, "data/selbal_HIV.csv")
readr::write_csv(as_tibble(selbal::sCD14), "data/selbal_sCD14.csv")
