
dir_init("./figures")

if (file.exists("raw_data.csv")) file.remove("raw_data.csv")
if (file.exists("prepped_data.csv")) file.remove("prepped_data.csv")

tic.clearlog()
