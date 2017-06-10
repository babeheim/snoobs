
rm(list=ls())

source('./code/project_functions.r')

my.data <- read.csv('./inputs/raw_data.csv', stringsAsFactors=FALSE)

prepped.data <- BDICE.sim.decomp.prepper(my.data, phenotype="snoob")

dir_init('./output')

write.csv(prepped.data, './output/prepped_data.csv', row.names=FALSE)

frq.table <- decomposer(prepped.data)

write.csv(frq.table, './output/frequency_table.csv', row.names=FALSE)

print("data decomposed")