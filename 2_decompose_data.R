
my.data <- read.csv("raw_data.csv", stringsAsFactors = FALSE)

prepped.data <- BDICE.sim.decomp.prepper(my.data, phenotype = "snoob")

write.csv(prepped.data, "prepped_data.csv", row.names = FALSE)

frq.table <- decomposer(prepped.data)

write.csv(frq.table, "figures/frequency_table.csv", row.names = FALSE)
