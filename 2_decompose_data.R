
cens <- read.csv("raw_data.csv")

deco <- BDICE.sim.decomp.prepper(cens, phenotype = "snoob")

# how do these tables differ? raw includes babies who died before their first census
# really that should be excluded since we dont know about them!
# deco does not include variables like age or sex!

deco <- as.data.frame(deco)

cens$key <- paste(cens$census, cens$id)
deco$key <- paste(deco$census, deco$id)

deco$age_days <- cens$age[match(deco$key, cens$key)]
deco$age <- round(deco$age_days / 365, 2)
deco$male <- cens$male[match(deco$key, cens$key)]
deco$snoob <- deco$phi

deco <- select(deco, -key)

write.csv(deco, "analysis_data.csv", row.names = FALSE)

frq.table <- decomposer(deco)
write.csv(frq.table, "figures/frequency_table.csv", row.names = FALSE)
