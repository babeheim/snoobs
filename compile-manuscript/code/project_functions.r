save_temp <- FALSE

my_seed <- 1701

dir_init <- function(path, verbose=FALSE){
    if(substr(path, 1, 2)!='./') stop('path argument must be formatted
        with "./" at beginning')
    contents <- dir(path, recursive=TRUE)
    if(verbose){
        if(length(contents)==0) print(paste('folder ', path, ' created.', sep=""))
        if(length(contents)>0) print(paste('folder ', path, ' wiped of ', length(contents), ' files/folders.', sep=""))
    }
    if(dir.exists(path)) unlink(path, recursive=TRUE)
    dir.create(path)
}

module_init <- function(path, verbose=FALSE){

    code_path <- paste(path, '/code', sep='')
    inputs_path <- paste(path, '/inputs', sep='')
    dir_init(path, verbose=verbose)
    dir_init(code_path, verbose=FALSE)
    dir_init(inputs_path, verbose=FALSE)

}

haversine <- function(lat1, lat2, long1, long2){

    R <- 6371000

    phi1 <- lat1*(2*pi)/360 # lat1 in radians
    phi2 <- lat2*(2*pi)/360  #lat2 in radians
    lambda1 <- long1*(2*pi)/360
    lambda2 <- long2*(2*pi)/360

    delta.phi <- phi2 - phi1  # lat2 - lat1 in radians
    delta.lambda <- lambda2 - lambda1  # in radians

    a <- sin(delta.phi/2)*sin(delta.phi/2) + cos(phi1)*cos(phi2)*sin(delta.lambda/2)*sin(delta.lambda/2)
    c <- 2*atan2(sqrt(a), sqrt(1-a))
    distance <- R*c
    return(distance)

}

modal <- function(data){
    mode <- names(sort(table(data),decreasing=T))[1]
    options(warn=-1)
    if(!is.na(as.numeric(mode))){
        mode <- as.numeric(mode)
    }
    options(warn=0)
    return(mode)
}

pid_maker <- function(n, reserved=""){
    my.let <- LETTERS
    my.num <- 0:9
    output <- replicate(n, paste(sample(c(my.let, my.num), 4, replace=TRUE), collapse=""))
    while(any(duplicated(output) | output %in% reserved)){
        output <- output[-which(output %in% reserved | duplicated(output))]
        remaining <- n-length(output)
        output <- c(output, replicate(remaining, paste(sample(c(my.let, my.num), 4, replace=TRUE), collapse="")))
    }
    return(output)
}

library(RColorBrewer)

kin_retriever <- function(pid, pop_reg){

    my_dad <- pop_reg$f_pid[pid]

    my_mom <- pop_reg$m_pid[pid]

    my_sons <- pop_reg$pid[which(pop_reg$f_pid==pid)]

    my_mates <- sort(unique(pop_reg$m_pid[which(pop_reg$f_pid==pid)]))

    my_mates_dads <- sort(unique(pop_reg$f_pid[which(pop_reg$pid %in% my_mates)]))

    my_mates_bros <- which(pop_reg$f_pid %in% my_mates_dads & pop_reg$male==1)
    
    my_sibs <- pop_reg$pid[which((pop_reg$f_pid==my_dad | pop_reg$m_pid==my_mom) & pop_reg$pid!=pid)]
    
    my_nephews <- pop_reg$pid[which(pop_reg$male==1 & (pop_reg$f_pid %in% my_sibs | pop_reg$m_pid %in% my_sibs))]
    
    my_grandparents <- sort(unique(c(pop_reg$ff_pid[pid], pop_reg$fm_pid[pid], pop_reg$mf_pid[pid], pop_reg$mm_pid[pid])))

    my_uncles <- integer(0)
    
    my_aunts <- integer(0)
    
    my_cousins <- integer(0)
    
    if(length(my_grandparents)>0){

        my_uncles <- pop_reg$pid[which(pop_reg$male==1 & (pop_reg$f_pid %in% my_grandparents | pop_reg$m_pid %in% my_grandparents) & pop_reg$pid!=my_dad)]

        my_aunts <- pop_reg$pid[which(pop_reg$male==0 & (pop_reg$f_pid %in% my_grandparents | pop_reg$m_pid %in% my_grandparents) & pop_reg$pid!=my_mom)]
        
        my_cousins <- pop_reg$pid[which((pop_reg$f_pid %in% my_uncles | pop_reg$m_pid %in% my_aunts) & pop_reg$pid!=pid)]

    }

    output <- na.omit(c(my_dad, my_mom, my_sons, my_mates, my_mates_dads, my_mates_bros, my_sibs, my_nephews, my_grandparents, my_uncles, my_aunts, my_cousins))

    output <- as.numeric(output)

    return(output)

}



neighbor_finder <- function(pid, pop_reg, house_reg){

    # for each pid, query their immediate neighbor network, including other men in their own house

    my_household <- pop_reg$household[pid]
    my_village <- pop_reg$village[pid]

    my_lat <- house_reg$lat[my_household]
    my_long <- house_reg$long[my_household]

    village_households <- which(house_reg$village==my_village)

    village_house_distances <- haversine(my_lat, hreg$lat[village_households], my_long, hreg$long[village_households])/1000
    nearby_houses <- village_households[which(village_house_distances<50)]

    my_neighbors <-  which(pop_reg$household %in% nearby_houses)
    my_neighbors <- setdiff(my_neighbors, pid)

    output <- as.numeric(na.omit(my_neighbors))

    return(output)

}

# diagnostic unit test for neighbor_finder
# test <- sample(active_rows,1)
# my_house <- preg$household[test] 
# my_village <-preg$village[test]
# plot(hreg$long[hreg$village==my_village], hreg$lat[hreg$village==my_village])
# points(hreg$long[my_house], hreg$lat[my_house], pch=20)
# my_neighbors <- neighbor_finder(test, preg, hreg)
# (my_neighbors)
# my_neighbors_houses <- preg$household[my_neighbors] 
# points(hreg$long[my_neighbors_houses], hreg$lat[my_neighbors_houses], pch=20, col='red')


spanishizer <- function(popreg){

    # everyone follows three trajectories 
    # rapid to fluency (rare)
    # minmal (the default)
    # market spanish (common)

    # what trajectory you are on is determined at birth, then you just follow
    # the trajectory throughout your life

    # spanish is a quantitative trait we treat as binary...?

    # if parents speak spanish
    # distance to SB
    # distance to school
    # if parents have a pequi

    # need to specify transition probs here...

}


pequi_get <- function(popreg, housereg, villagereg){

    logit <- function (x) log(x) - log(1 - x)
    logistic <- function (x) 1/(1+exp(-x))

    can <- which(is.na(popreg$dod) & popreg$male==1 & popreg$age >= 16*365 & popreg$pequi==0)

    dist_SB <- villagereg$dist_SB[popreg$village[can]]/1000/10
    kin_has <- rep(0, length(can))
    neighbor_has <- rep(0, length(can))

    for(i in 1:length(can)){

        my_kin <- kin_retriever(i, popreg)
        kin_has[i] <- any(popreg$pequi[my_kin]==1, na.rm=TRUE)

        my_neighbors <- neighbor_finder(i, popreg, housereg)
        neighbor_has[i] <- any(popreg$pequi[my_neighbors]==1, na.rm=TRUE)

    }

    baseline <- 0.0001 # [0,1]
    alpha <- log(baseline/(1-baseline))

    beta_neighbor <- 3
    beta_kin <- -1
    beta_dist <- 0.5

    logit_pr_pequi_get <-
        alpha + 
        beta_neighbor * neighbor_has + 
        beta_kin * kin_has + 
        beta_dist * dist_SB

    got_pequi <- rbinom(length(can), 1, logistic(logit_pr_pequi_get))

    output <- can[which(got_pequi==1)]

    # # diagnostics - is our children learning?
    
    # if(runif(1) < 0.1){
       #  # tar <- which(got_pequi==1)
       #  tar <- 1:length(can)
       #  o <- cbind(can[tar], dist_SB[tar], kin_has[tar], neighbor_has[tar], logit_pr_pequi_get[tar], got_pequi[tar])
       #  colnames(o) <- c('pid', 'dist_SB', 'kin_has', 'neighbor_has', 'logit_pr', 'pequi')
       #  o <- as.data.frame(o)
       #  # lm(logit_pr ~ dist_SB + kin_has + neighbor_has, data=o)
       #  write.csv(o, './temp/pequi_snoop.csv', row.names=FALSE)
    # }

    return(output)

}

daily_immigration <- function(crude_immigration_rate_f, data){
    daily_immigration_rate_per_1000 <- (crude_immigration_rate_f/1000)/365
    pop_size <- length(active_rows)
    mean_daily_n_immigrants <- daily_immigration_rate_per_1000*pop_size 
    todays_n_immigrants <- rpois(1, mean_daily_n_immigrants)        
    unoccupied_rows <- sum(is.na(data[,1]))
    return(todays_n_immigrants)
}

conception <- function(asf, register){

    active_rows <- which(is.na(register$dod))

    male <- register$male[active_rows]
    not_pregnant <- is.na(register$counter[active_rows])
    active_ages <- register$age[active_rows]
    fecund_women_rows <- active_rows[which(male == 0 & not_pregnant & 
        active_ages >= 15*365 & active_ages < 50*365)]

    snoob <- register[fecund_women_rows,'snoob']
    
    asfr_years <- asf$asfr/c(diff(asf$age), 1)

    baseline <- (asfr_years[findInterval(register$age[fecund_women_rows]/365, 
        asf$age)]/1000)/365
    alpha = log(baseline/(1- baseline))
    
    daily_pr_conception <- (1-(1/(1+exp(alpha + 
        fertility_bias_snoob*register$snoob[fecund_women_rows]))))
    got_pregnant <- rbinom(length(fecund_women_rows), 1, daily_pr_conception)
    return(fecund_women_rows[which(got_pregnant==1)])

}

grim_reaper <- function(asm, register){

    active_rows <- which(is.na(register$dod))

    asmr_years <- asm$asmr/c(diff(asm$age), 1)

    baseline <- asmr_years[findInterval(register$age[active_rows]/365, 
        asm$age)]/1000/365
    alpha = log(baseline/(1 - baseline))

    daily_pr_death <- (1-(1/(1+exp(alpha + 
        mortality_bias_male*register$male[active_rows] + 
        mortality_bias_snoob*register$snoob[active_rows]))))
    died <- rbinom(length(active_rows), 1, daily_pr_death)
    return(active_rows[which(died==1)])
    
}

baby_maker <- function(mom_id, dad_id, register){
    new <- register[1,]
    new[!is.na(new)] <- NA
    new$pid <- nrow(register)+1
    new$age <- 0
    new$dob <- day # scoping problem___
    new$male <- rbinom(1,1,prob=0.5)
    new$m_pid <- mom_id
    new$f_pid <- dad_id

    new$ff_pid <- register$f_pid[dad_id]
    new$fm_pid <- register$m_pid[dad_id]
    new$mf_pid <- register$f_pid[mom_id]
    new$mm_pid <- register$m_pid[mom_id]

    new$household <- register$household[mom_id]
    new$village <- register$village[mom_id]
    new$pequi <- 0

    # snoob inheritance
    mom_snoob <- register$snoob[mom_id]
    dad_snoob <- register$snoob[dad_id]
    midparent_snoob <- mean(c(mom_snoob, dad_snoob))  
    # GOD DAMN that was a hard bug to solve___
    if(midparent_snoob == 1) midparent_snoob <- 0.9999
    if(midparent_snoob == 0) midparent_snoob <- 0.0001
    baseline_trans_bias <- midparent_snoob
    alpha <- log(baseline_trans_bias/(1-baseline_trans_bias))
    prob_kid_is_snoob <- exp(alpha + transmission_bias_snoob)/
        (1+ exp(alpha + transmission_bias_snoob))
    new$snoob <- rbinom(1,1,prob=prob_kid_is_snoob)
    
    # english inheritance
    mom_english <- register$language[mom_id]
    dad_english <- register$language[dad_id]
    parent_frq_english <- mean(c(mom_english, dad_english))
    parent_alpha <- 0.99
    pop_frq_english <- mean(register$language[active_rows])
    pr_acquire_english <- parent_alpha*parent_frq_english + 
        (1-parent_alpha)*pop_frq_english
    new$language <- rbinom(1, 1, pr_acquire_english)
    
    new$h_gene <- sample( c(register$h_gene[mom_id], 
        register$h_gene[mom_id]) , 1)
    new$height <- height_fun(age=new$age, h_gene=new$h_gene) 
    
    register <- rbind(register, new)
}

mate_finder <- function(mom_id, data){
    
    current_mate <- NA
    
    moms_snoob <- data[mom_id, 'snoob']
    
    last_mate <- data[mom_id, 'mate']
    if(!is.na(last_mate)){
        if(is.na(data$dod[last_mate]) & data[last_mate, 'age'] > 15*365 & 
            data[last_mate, 'age'] < 60*365){
            current_mate <- last_mate
        }
    }
    # if the last mate exists, is alive and within the age range, he's your man
    
    if(is.na(current_mate)){    
        # 'available' means within age range, alive, male
        # male_mating_bias <- 1
        
        available_men <- active_rows[which( data[active_rows, 'male']==1 & 
            data[active_rows, 'age'] >= 15*365 & 
            data[active_rows, 'age'] < 60*365 )]
                    
        if(length(available_men)>0){
            current_mate <- available_men[1]
            if(length(available_men) > 1){
                available_men_traits <- data[available_men,'snoob']
                n_men <- length(available_men)
            
                prob_choose_snoob <- mate_similarity_bias_snoob*moms_snoob + 
                    (1-moms_snoob)*(1-mate_similarity_bias_snoob)
                prob_choose_nonsnoob <- 1 - prob_choose_snoob
    
                prob_is_chosen <- rep(0.5, n_men)
                prob_is_chosen[available_men_traits==1] <- prob_choose_snoob
                prob_is_chosen[available_men_traits==0] <- prob_choose_nonsnoob
                                
                current_mate <- sample(available_men, 1, prob=prob_is_chosen)
            }
        } else {
            current_mate <- NA
        }
    }
    current_mate

}

height_fun <- function(age, h_gene){
    if(any(age>100)) age <- age/365
    heights <- h_gene*exp(-3 + .2*age)/(1 + exp(-3+.2*age))
    heights
}

group_maker <- function(n, reserved='', seed=NA){
    my_let <- LETTERS 
    my_num <- 0:9 
    if(is.na(seed) | !is.numeric(seed)) set.seed(as.numeric(as.POSIXlt(Sys.time())))
    if(!is.na(seed) & is.numeric(seed)) set.seed(seed)
    output <- replicate(n, paste(sample(c(my_let, my_num), 4, replace=TRUE), 
        collapse=''))
    while(any(duplicated(output) | output %in% reserved)){
        output <- output[-which(output %in% reserved | duplicated(output))]
        remaining <- n-length(output)
        output <- c(output, replicate(remaining, paste(sample(c(my_let, my_num), 4, 
            replace=TRUE), collapse='')))
    }
    output
}

