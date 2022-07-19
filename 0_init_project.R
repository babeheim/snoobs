
dir_init("./figures")

if (file.exists("censuses.csv")) file.remove("censuses.csv")
if (file.exists("analysis_data.csv")) file.remove("analysis_data.csv")

tic.clearlog()
