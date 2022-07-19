
[x] flatten repo
[x] add docstrings
[x] rename variables to use dollar signs
[x] remove all commented-out code
[x] add stopifnot checks
[x] refactor analysis code
[x] remove unnecessary features of the simulation
[x] remove hard-coded coefficients from figure code
[x] replace calls of `sample.naive.posterior`
[x] rename 'state' in the simulator to 'is_alive'
[x] remove reference of a 'state' variable in phisim
[x] remove intercensus ghosts from the raw data before prepping
[x] index census from 1 to n_censuses
[x] reassign codes to be more natural
[x] update all the census and prepped data variable names
[x] re-write decomposer so it doesnt need the dead persons phenotype for the census where they are already recorded as dead
[x] re-write this to work from 1:(n_censuses - 1) as in the worked table 2 example
[x] refactor `prep_census_data` & bring into the main script
[x] refactor `decompose_prepped_census` & bring into the main script
[x] remove 'grouping' aspect of `decompose_prepped_census`
[x] rewrite the census-recoding to ignore people born and died intercensus
[x] relabel the censuses to start at 1 not 0
[x] confrim the `n.kids` in the prepped data is the kids you WILL have next intercensus
[x] rewrite analysis script to use the new nomenclature and variable definitions

[-] update the table 2 work to incorporate my new insights into the right way to store info, thinking of pairs of entries

[-] think up better name than `analysis_data.csv`
[-] refactor fertility and mortality figures
[-] fix the fit on the individual change figure
[-] flag all warnings
[-] output supplementary tables as tex files
[-] why does my analysis figure code assume that phi will be near 0.5 at some point?
