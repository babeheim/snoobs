
cens <- read.csv("raw_data.csv")

prep <- prep_census_data(cens, phenotype = "snoob")

# how do these tables differ? raw includes babies who died before their first census
# really that should be excluded since we dont know about them!
# prep does not include variables like age or sex!

prep <- as.data.frame(prep)

cens$key <- paste(cens$census, cens$id)
prep$key <- paste(prep$census, prep$id)

prep$age_days <- cens$age[match(prep$key, cens$key)]
prep$age <- round(prep$age_days / 365, 2)
prep$male <- cens$male[match(prep$key, cens$key)]
prep$snoob <- prep$phi

prep <- select(prep, -key)

write.csv(prep, "analysis_data.csv", row.names = FALSE)

frq.table <- decompose_prepped_census(prep)
write.csv(frq.table, "figures/frequency_table.csv", row.names = FALSE)
