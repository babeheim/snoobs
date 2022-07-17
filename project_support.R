
library(RColorBrewer)
library(bbmle)
library(rethinking)
library(tictoc)

machine_name <- "thinkpad x360"
project_seed <- 1912
set.seed(project_seed)

source("R/misc_functions.R")
source("R/simulation_functions.R")
source("R/decomposition_functions.R")
source("R/analysis_functions.R")
