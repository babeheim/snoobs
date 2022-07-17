
rm(list = ls())

source("./project_support.R")

dir_init("./figures")

if (file.exists("raw_data.csv")) file.remove("raw_data.csv")
if (file.exists("prepped_data.csv")) file.remove("prepped_data.csv")

tic.clearlog()

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
  task = task,
  time_min = time_min
)

write.csv(report, file.path("figures/timing-report.csv"), row.names = FALSE)
