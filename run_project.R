
rm(list = ls())

source("project_support.R")

source("0_init_project.R")

##############

tic("run snoobsim decomposition analysis")

tic("simulate data")
source("1_simulate_data.R")
toc(log = TRUE)

tic("decompose simulation data")
source("2_decompose_data.R")
toc(log = TRUE)

tic("analyze data")
source("3_analyze_data.R")
toc(log = TRUE)

toc(log = TRUE)

###########

tic.log(format = TRUE)
msg_log <- unlist(tic.log())

task <- msg_log
task <- gsub(":.*$", "", task)

time_min <- msg_log
time_min <- gsub("^.*: ", "", time_min)
time_min <- gsub(" sec elapsed", "", time_min)
time_min <- round(as.numeric(time_min)/60, 2)

report <- data.frame(
  project_seed = project_seed,
  machine = machine_name,
  n_sim_years = n_sim_years,
  n_ind_init = n_ind_init,
  census_interval_years = census_interval_years,
  frq_snoob_init = frq_snoob_init,
  task = task,
  time_min = time_min
)

write.csv(report, file.path("figures/timing-report.csv"), row.names = FALSE)
