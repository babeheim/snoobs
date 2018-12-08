
rm(list = ls())
start_time <- Sys.time()
source("./project_support.r")
tic.clearlog()

##############

tic("simulate data")
setwd("./1_simulate_data")
source("./simulate_data.r")
setwd("..")
toc(log = TRUE)

##############

tic("decompose simulation data")
dir_init("./2_decompose_data/inputs")
file.copy("./1_simulate_data/output/raw_data.csv", "./2_decompose_data/inputs")
setwd("./2_decompose_data")
source("./decompose_data.r")
setwd("..")
toc(log = TRUE)

##############

tic("analyze data")
dir_init("./3_analyze_data/inputs")
file.copy("./1_simulate_data/output/raw_data.csv", "./3_analyze_data/inputs")
file.copy("./2_decompose_data/output/prepped_data.csv", "./3_analyze_data/inputs")
setwd("./3_analyze_data")
source("./analyze_data.r")
setwd("..")
toc(log = TRUE)

##############

tic("make manuscript")
dir_init("./4_make_manuscript/inputs")
setwd("./4_make_manuscript")
source("./make_manuscript.r")
setwd("..")
toc(log = TRUE)

##############

dir_init("./output")
if (!exists("start_time")) start_time <- "unknown"
write_log(title = "project: SnoobSim decomposition analysis",
  path = "./output/log.txt", start_time = start_time)
