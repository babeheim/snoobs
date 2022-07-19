
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
[-] think up better id than `analysis_data.csv`
[-] refactor `prep_census_data`
[-] refactor `decompose_prepped_census`
[-] remove 'grouping' aspect of `decompose_prepped_census`
[-] edge case: people born but parents were not present in the population last census are immigrants, not borns
[-] refactor fertility and mortality figures
[-] rewrite the census-recoding to ignore people born and died intercensus
[-] relabel the censuses to start at 1 not 0
[-] fix the fit on the individual change figure
[-] confrim the `n.kids` in the prepped data is the kids you WILL have next intercensus
[-] flag all warnings
[-] output supplementary tables as tex files
[-] why does my analysis figure code assume that phi will be near 0.5 at some point?
