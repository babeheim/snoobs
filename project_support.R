
library(dplyr)
library(RColorBrewer)
library(bbmle)
library(rethinking)
library(tictoc)

n_sim_years <- 300 # publication value 300 years
n_ind_init <- 1000 # publication value 1000 people
census_interval_years <- 5 # publication value 5 years
frq_snoob_init <- 0.2

project_seed <- 1912
set.seed(project_seed)
machine_name <- "thinkpad x360"

source("R/misc_functions.R")
source("R/simulation_functions.R")
source("R/analysis_functions.R")
