
dir_init("./figures")

if (file.exists("raw_data.csv")) file.remove("raw_data.csv")
if (file.exists("analysis_data.csv")) file.remove("analysis_data.csv")

tic.clearlog()
