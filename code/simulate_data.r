
####### Simulation Control Panel #######

source('./code/project_functions.r')

# Basic Simulator Values
census.census.length <- 365*5 # number of days until we do another census
diagnostic = FALSE
days.max = 365*300 # stop after this time
running.pop.graphs = TRUE
event.logging = FALSE
reproductive.min.age <- 14 # in years
ini.snoob.frq <- 0.2
ini.pop.size <- 1000

# Vital Population Rates
age.specific.fertility <- c(40.47, 102.50, 115.92, 96.10, 46.36, 9.12, 0.58) # expected births per person within five-year age classes from 15 to 50
disable.births <- FALSE # for diagnostic purposes
age.specific.mortality <- c(6.54, 0.29, 0.14, 0.18, 0.64, 0.91, 0.90, 1.06, 1.53, 2.31, 3.41, 4.93, 7.42, 11.50, 17.80, 27.71, 43.50, 69.58, 110.56, 174.77, 276.56, 438.92) # expected deaths per person within five-year age classes from 0 to 100+
disable.deaths <- FALSE  # for diagnostic purposes
crude.immigration.rate <- 0 # expected number of immigrants per 1000 people per year     
crude.emigration.rate <-  0 # expected number of emigrants per 1000 people per year
individual.change.rate <- 0.00006  # fraction of individuals who update each day
social.sample.size <- 100  # number of individuals sampled when updating
mortality.bias.male <- 0.3 # logistic coefficient: positive values mean men die more often

# Mechanisms of Evoluton
mate.similarity.bias.snoob <- 0.9 # from 0 to 1; if 0.5, no assortment by snoob status.  If 1, perfect assortment along snoob status.  If 0, perfect reverse assortment along snoob status.
fertility.bias.snoob <- 0.8 # logistic coefficient; positive values mean snoobs have higher RS than nonsnoobs, all else being equal
immigration.bias.snoob <- 0.0  # logistic coefficient; positive values mean immigrants are more likely to be snoobs than would be expected for the current frequency of snoobs in the population.  When zero, immigrants contribute trait at same rate as current frq of snoobs.
emigration.bias.snoob <- 0.0  # logistic coefficient: positive values mean snoobs emigrate more often.  
mortality.bias.snoob <- -0.5 # logistic coefficient: positive values mean snoobs die more often.
transmission.bias.snoob <- 0.0   # logistic coefficient: if 0, kid will take snoob=1 with probability equal to midparent value.  If >0, kid will be more likely to take snoob=1 than would be expected from the midparent value.  If <0, kid will be LESS likely to take snoob=1 than would be expected from midparent values....to see real effects, this term has to be LARGE - magnitude 10 or so.
conformity.bias.snoob <- 1.3 # exponential coefficient: unbiased transmission if equals 1....values > 1 mean biased copying towards majority's trait, values < 1 mean biased copying against majority's trait
snoob.less.risky.factor <- 0.9 # varies between 0 and 1; 1 indicates snoobs are just as risky as non-snoobs, while 0 indicates they aren't risky at all

# Diagnostic on parameters; if there's a problem, parameters will not be okay.  
parameters.okay <- TRUE
if(mate.similarity.bias.snoob > 1 | mate.similarity.bias.snoob < 0) parameters.okay <- FALSE
parameters.okay


####### Simulator Initialization #######

# SnoobSim State Codes:
# 0 - Dead 
# 1 - Active 
# 2 - Emigrant

reg.max <- ini.pop.size*2  # this isn't really important anymore...the matrix will grow

day <- 1
x.axis.day <- 1

standard.traits <- c("state", "mom", "dad", "age", "age.cat", "last.im", "last.em", "died", "male",  "mate", "counter")
special.traits <- c("h.gene", "height", "snoob", "language", "risk")

census.n.opening <- ini.pop.size
census.n.births <- 0
census.n.deaths <- 0
census.n.immigrants <- 0
census.n.emigrants <- 0
census.n.closing <- 0

census.table <- c(census.n.opening, census.n.births, census.n.deaths, census.n.immigrants, census.n.emigrants, census.n.closing)
names(census.table) <- c("n.start", "births", "deaths", "inflow", "outflow", "n.close")

state <- c(rep(1, ini.pop.size), rep(NA, reg.max-ini.pop.size))
mom <- c(rep(0, ini.pop.size), rep(NA, reg.max-ini.pop.size))

pop.reg <- matrix(NA,ncol=length(c(standard.traits, special.traits)), nrow=reg.max)
pop.reg <- as.data.frame(pop.reg)
colnames(pop.reg) <- c(standard.traits, special.traits)

pop.reg[1:ini.pop.size, "state"] <- 1
pop.reg[1:ini.pop.size, "mom"] <- 0
pop.reg[1:ini.pop.size, "dad"] <- 0


ages <- 0:99
counts <- c(57, 74, 68, 77, 65, 75, 79, 83, 66, 60, 80, 72, 59, 87, 78, 73, 58, 73, 84, 67, 83, 76, 71, 74, 72, 62, 91, 64, 63, 76, 80, 79, 71, 73, 76, 72, 70, 82, 80, 70, 82, 72, 78, 75, 75, 83, 73, 70, 79, 71, 77, 81, 80, 71, 64, 70, 62, 65, 73, 76, 70, 73, 84, 74 ,64 ,63 ,76 ,72 ,55 ,81 ,78 ,74 ,68 ,69 ,71 ,63 ,56, 50, 55, 12, 41, 38, 56, 30, 27, 22, 38, 30,  8, 22, 11 ,12, 12, 10,  4,  9,  8,  3,  2,  4)

pop.reg[1:ini.pop.size, "age"] <- 365*sample(ages, ini.pop.size, replace=T, prob=counts)+replicate(ini.pop.size, sample(0:364,1))
# pop.reg[1:ini.pop.size, "age"] <- 365*sample(1:60, ini.pop.size, replace=T, prob=dexp(1:60, rate=.05))
pop.reg[1:ini.pop.size, "age.cat"] <- age.binner(pop.reg[1:ini.pop.size, "age"])
pop.reg[1:ini.pop.size, "last.im"] <- pop.reg[1:ini.pop.size, "age"]
pop.reg[1:ini.pop.size, "male"] <- rep(c(1,0), ini.pop.size/2)

# trait initialization
pop.reg[1:ini.pop.size, "h.gene"] <- sample( 80:120 , ini.pop.size, replace=T)
pop.reg[1:ini.pop.size, "height"] <- height.fun(pop.reg[1:ini.pop.size, "age"], pop.reg[1:ini.pop.size, "h.gene"])
pop.reg[1:ini.pop.size, "snoob"] <- rbinom(ini.pop.size, 1, p=ini.snoob.frq)
pop.reg[1:ini.pop.size, "language"] <- pop.reg[1:ini.pop.size, "snoob"]
# uncorrelate.switch <- sample(1:ini.pop.size, 0.15*ini.pop.size) # If I did this right, snoobs and language should correlate very nicely
# pop.reg[uncorrelate.switch, "language"] <- 1-pop.reg[uncorrelate.switch, "language"]
pop.reg[1:ini.pop.size, "risk"] <- rep(0, ini.pop.size)

na.reg <- matrix(NA,ncol=length(c(standard.traits, special.traits)), nrow=reg.max*.3)

census.death.list <- integer(0)
census.emigration.list <- integer(0)

event.log <- character(0)

avg.kiddy <- integer(0)

n <- ini.pop.size
n.snoob <- sum(pop.reg[1:ini.pop.size, "snoob"])
n.eng <- sum(pop.reg[1:ini.pop.size, "language"])
# plot(day, n, pch=20, ylim=c(0, reg.max), xlim=c(0, days.max))

if(disable.births == T) age.specific.fertility <- rep(0, 7)
if(disable.deaths == T) age.specific.mortality <- rep(0, 22)
active.rows <- which(pop.reg[,"state"]==1)
pop.reg[active.rows, "risk"] <- risk.setter(pop.reg)

census.end.census <- cbind(rep(0, length(active.rows)), active.rows, pop.reg[active.rows,])
colnames(census.end.census)[1] <- "census"
colnames(census.end.census)[2] <- "id"

names(age.specific.fertility) <- c(15, 20, 25, 30, 35, 40, 45)
names(age.specific.mortality) <- c(0, 1, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 100)

should.simulate.day <- day < days.max & sum(is.na(pop.reg[,1])) != 0




diagnostics <- F

if(diagnostics==TRUE){
daily.age.counts <- table(floor(pop.reg[active.rows, "age"]/365))
M.at.risk <- rep(0, 130)
names(M.at.risk) <- 0:129
M.events <- rep(0, 130)
names(M.events) <- 0:129

snoob.daily.age.counts <- daily.age.counts
M.snoob.at.risk <- rep(0, 130)
names(M.snoob.at.risk) <- 0:129
M.snoob.events <- rep(0, 130)
names(M.snoob.events) <- 0:129

nonsnoob.daily.age.counts <- daily.age.counts
M.nonsnoob.at.risk <- rep(0, 130)
names(M.nonsnoob.at.risk) <- 0:129
M.nonsnoob.events <- rep(0, 130)
names(M.nonsnoob.events) <- 0:129

snooblish.daily.age.counts <- daily.age.counts
M.snooblish.at.risk <- rep(0, 130)
names(M.snooblish.at.risk) <- 0:129
M.snooblish.events <- rep(0, 130)
names(M.snooblish.events) <- 0:129

nonsnooblish.daily.age.counts <- daily.age.counts
M.nonsnooblish.at.risk <- rep(0, 130)
names(M.nonsnooblish.at.risk) <- 0:129
M.nonsnooblish.events <- rep(0, 130)
names(M.nonsnooblish.events) <- 0:129

nonsnoob.snooblish.daily.age.counts <- daily.age.counts
M.nonsnoob.snooblish.at.risk <- rep(0, 130)
names(M.nonsnoob.snooblish.at.risk) <- 0:129
M.nonsnoob.snooblish.events <- rep(0, 130)
names(M.nonsnoob.snooblish.events) <- 0:129

nonsnoob.nonsnooblish.daily.age.counts <- daily.age.counts
M.nonsnoob.nonsnooblish.at.risk <- rep(0, 130)
names(M.nonsnoob.nonsnooblish.at.risk) <- 0:129
M.nonsnoob.nonsnooblish.events <- rep(0, 130)
names(M.nonsnoob.nonsnooblish.events) <- 0:129

snoob.snooblish.daily.age.counts <- daily.age.counts
M.snoob.snooblish.at.risk <- rep(0, 130)
names(M.snoob.snooblish.at.risk) <- 0:129
M.snoob.snooblish.events <- rep(0, 130)
names(M.snoob.snooblish.events) <- 0:129

snoob.nonsnooblish.daily.age.counts <- daily.age.counts
M.snoob.nonsnooblish.at.risk <- rep(0, 130)
names(M.snoob.nonsnooblish.at.risk) <- 0:129
M.snoob.nonsnooblish.events <- rep(0, 130)
names(M.snoob.nonsnooblish.events) <- 0:129

}



####### Simulator Daily Loop #######

while( should.simulate.day ){
        
    if(day %% census.census.length == 1) census.n.opening <- sum(pop.reg[,"state"]==1, na.rm=T)

    active.rows <- which(pop.reg[,"state"]==1)
    
    # set risk and age category for the day
    pop.reg[active.rows, "risk"] <- risk.setter(pop.reg)
    pop.reg[active.rows, "age.cat"] <- age.binner(pop.reg[active.rows, "age"])
    
    # immigration
    n.immigrants <- daily.immigration(crude.immigration.rate, pop.reg)
    census.n.immigrants <- census.n.immigrants + n.immigrants
    if(n.immigrants > 0){
        for(i in 1:n.immigrants){
            im.surname <- min(which(is.na(pop.reg[,1])))
            pop.reg[im.surname, ] <- immigrant.maker(data=pop.reg)
            if(is.na(pop.reg[im.surname, "snoob"])) stop("WTF?")
            if(event.logging==TRUE){
                event.log <- c(event.log, paste("Day ", day, ": ", ifelse(pop.reg[im.surname, "male"]==1, "Male ", "Female "), im.surname, " has entered the population.", sep=""))
            }
        }
    }
    
    active.rows <- which(pop.reg[,"state"]==1)
    
    # births
    giving.birth <- active.rows[!is.na(pop.reg[active.rows,"counter"]) & pop.reg[active.rows,"counter"]==0]
    if(length(giving.birth)>0){
        for(i in 1:length(giving.birth)){
            mom <- giving.birth[i]
            dad <- pop.reg[giving.birth[i], "mate"]
            baby.surname <- min(which(is.na(pop.reg[,1])))
            pop.reg[baby.surname, ] <- baby.maker(mom, dad, pop.reg)
            if(event.logging==TRUE){
                event.log <- c(event.log, paste("Day ", day, ": ", ifelse(pop.reg[baby.surname, "male"]==1, "Male ", "Female "), baby.surname, " has been born to female ", pop.reg[baby.surname, "mom"] , " and male ", pop.reg[baby.surname, "dad"], ".", sep=""))
            }
            census.n.births <- census.n.births + 1    
            pop.reg[mom, "counter"] <- NA
        }
    }
    
    if(diagnostics==TRUE){
        # at risk is any woman in the age range
        # event is any woman who just gave birth
    }
    
    active.rows <- which(pop.reg[,"state"]==1)
    
    # mating
    male <- pop.reg[active.rows, "male"]
    not.pregnant <- is.na(pop.reg[active.rows, "counter"])
    active.age.cats <- pop.reg[active.rows, "age.cat"]
    fecund.women.rows <- active.rows[which(male == 0 & not.pregnant & active.age.cats >= 15 & active.age.cats < 50)]
    will.conceive <- conception(age.specific.fertility, pop.reg)
    female.maters <- fecund.women.rows[as.logical(will.conceive)]
    female.maters

    # diagnostics: make sure only women are in fecund.women.rows
    
    if(length(female.maters) > 0){
        for(i in 1:length(female.maters)){
            mom <- female.maters[i]
            dad <- mate.finder(mom, data=pop.reg)
            pop.reg[mom,"mate"] <- dad
            if(is.na(dad)){
                event.log <- c(event.log, paste("Day ", day, ": ", "Female ", mom, " tried to conceive but was unable to find a mate.", sep=""))
            }
            if(!is.na(dad)){
                pop.reg[dad,"mate"] <- mom
                pop.reg[mom, "counter"] <- round(rnorm(1,280,5))
                if(event.logging==TRUE){
                    event.log <- c(event.log, paste("Day ", day, ": ", "Female ", mom, " has become pregnant from male ", dad, ".", sep=""))
                }
            }        
        }
    }
    
    # mortality
    who.died <- grim.reaper(age.specific.mortality, data=pop.reg)
    census.n.deaths <- census.n.deaths + sum(who.died)    
    who.died <- as.logical(who.died)
        
    if(any(who.died)){
        census.death.list <- c(census.death.list, active.rows[who.died])
        pop.reg[active.rows[who.died], "state"] <- 0
        pop.reg[active.rows[who.died],"died"] <- pop.reg[active.rows[who.died], "age"]  # record their age at death
        if(event.logging==TRUE){
            death.list <- active.rows[who.died]
            for(i in 1:length(death.list)){
                event.log <- c(event.log, paste("Day ", day, ": ", ifelse(pop.reg[death.list[i], "male"]==1, "Male ", "Female "), death.list[i], " has died.", sep=""))
            }
        }
    }
    
    if(diagnostics==TRUE){
        daily.age.counts <- table(floor(pop.reg[active.rows, "age"]/365))
         M.at.risk[names(daily.age.counts)] <- M.at.risk[names(daily.age.counts)]+daily.age.counts
        if(any(who.died)){
            daily.event.age <- table(floor(pop.reg[active.rows[who.died], "age"]/365))
            M.events[names(daily.event.age)] <-  M.events[names(daily.event.age)]+daily.event.age
        }
        
        is.snoob <- as.logical(pop.reg[active.rows,"snoob"])
        
        snoob.daily.age.counts <- table(floor(pop.reg[active.rows[is.snoob], "age"]/365))
        M.snoob.at.risk[names(snoob.daily.age.counts)] <- M.snoob.at.risk[names(snoob.daily.age.counts)]+snoob.daily.age.counts
        if(any(is.snoob & who.died)){
            snoob.daily.event.age <- table(floor(pop.reg[active.rows[which(is.snoob & who.died)], "age"]/365))
            M.snoob.events[names(snoob.daily.event.age)] <- M.snoob.events[names(snoob.daily.event.age)]+snoob.daily.event.age
        }
        
        nonsnoob.daily.age.counts <- table(floor(pop.reg[active.rows[!is.snoob], "age"]/365))
        M.nonsnoob.at.risk[names(nonsnoob.daily.age.counts)] <- M.nonsnoob.at.risk[names(nonsnoob.daily.age.counts)]+nonsnoob.daily.age.counts
        if(any(!is.snoob & who.died)){
            nonsnoob.daily.event.age <- table(floor(pop.reg[active.rows[which(!is.snoob & who.died)], "age"]/365))
            M.nonsnoob.events[names(nonsnoob.daily.event.age)] <- M.nonsnoob.events[names(nonsnoob.daily.event.age)]+nonsnoob.daily.event.age
        }
        
        is.snooblish <- as.logical(pop.reg[active.rows,"language"])
        
        snooblish.daily.age.counts <- table(floor(pop.reg[active.rows[is.snooblish], "age"]/365))
        M.snooblish.at.risk[names(snooblish.daily.age.counts)] <- M.snooblish.at.risk[names(snooblish.daily.age.counts)]+snooblish.daily.age.counts
        if(any(is.snooblish & who.died)){
            snooblish.daily.event.age <- table(floor(pop.reg[active.rows[which(is.snooblish & who.died)], "age"]/365))
            M.snooblish.events[names(snooblish.daily.event.age)] <- M.snooblish.events[names(snooblish.daily.event.age)]+snooblish.daily.event.age
        }
        
        nonsnooblish.daily.age.counts <- table(floor(pop.reg[active.rows[!is.snooblish], "age"]/365))
        M.nonsnooblish.at.risk[names(nonsnooblish.daily.age.counts)] <- M.nonsnooblish.at.risk[names(nonsnooblish.daily.age.counts)]+nonsnooblish.daily.age.counts
        if(any(!is.snooblish & who.died)){
            nonsnooblish.daily.event.age <- table(floor(pop.reg[active.rows[which(!is.snooblish & who.died)], "age"]/365))
            M.nonsnooblish.events[names(nonsnooblish.daily.event.age)] <- M.nonsnooblish.events[names(nonsnooblish.daily.event.age)]+nonsnooblish.daily.event.age
        }
        
        snoob.snooblish.daily.age.counts <- table(floor(pop.reg[active.rows[is.snooblish & is.snoob], "age"]/365))
        M.snoob.snooblish.at.risk[names(snoob.snooblish.daily.age.counts)] <- M.snoob.snooblish.at.risk[names(snoob.snooblish.daily.age.counts)] + snoob.snooblish.daily.age.counts
        if(any(is.snooblish & is.snoob & who.died)){
            snoob.snooblish.daily.event.age <- table(floor(pop.reg[active.rows[which(is.snooblish & is.snoob & who.died)], "age"]/365))
            M.snoob.snooblish.events[names(snoob.snooblish.daily.event.age)] <- M.snoob.snooblish.events[names(snoob.snooblish.daily.event.age)]+snoob.snooblish.daily.event.age
        }
        
        snoob.nonsnooblish.daily.age.counts <- table(floor(pop.reg[active.rows[!is.snooblish & is.snoob], "age"]/365))
        M.snoob.nonsnooblish.at.risk[names(snoob.nonsnooblish.daily.age.counts)] <- M.snoob.nonsnooblish.at.risk[names(snoob.nonsnooblish.daily.age.counts)] + snoob.nonsnooblish.daily.age.counts
        if(any(!is.snooblish & is.snoob & who.died)){
            snoob.nonsnooblish.daily.event.age <- table(floor(pop.reg[active.rows[which(!is.snooblish & is.snoob & who.died)], "age"]/365))
            M.snoob.nonsnooblish.events[names(snoob.nonsnooblish.daily.event.age)] <- M.snoob.nonsnooblish.events[names(snoob.nonsnooblish.daily.event.age)]+snoob.nonsnooblish.daily.event.age
        }
        
        nonsnoob.nonsnooblish.daily.age.counts <- table(floor(pop.reg[active.rows[!is.snooblish & !is.snoob], "age"]/365))
        M.nonsnoob.nonsnooblish.at.risk[names(nonsnoob.nonsnooblish.daily.age.counts)] <- M.nonsnoob.nonsnooblish.at.risk[names(nonsnoob.nonsnooblish.daily.age.counts)] + nonsnoob.nonsnooblish.daily.age.counts
        if(any(!is.snooblish & !is.snoob & who.died)){
            nonsnoob.nonsnooblish.daily.event.age <- table(floor(pop.reg[active.rows[which(!is.snooblish & !is.snoob & who.died)], "age"]/365))
            M.nonsnoob.nonsnooblish.events[names(nonsnoob.nonsnooblish.daily.event.age)] <- M.nonsnoob.nonsnooblish.events[names(nonsnoob.nonsnooblish.daily.event.age)]+nonsnoob.nonsnooblish.daily.event.age
        }
        
        nonsnoob.snooblish.daily.age.counts <- table(floor(pop.reg[active.rows[is.snooblish & !is.snoob], "age"]/365))
        M.nonsnoob.snooblish.at.risk[names(nonsnoob.snooblish.daily.age.counts)] <- M.nonsnoob.snooblish.at.risk[names(nonsnoob.snooblish.daily.age.counts)] + nonsnoob.snooblish.daily.age.counts
        if(any(is.snooblish & !is.snoob & who.died)){
            nonsnoob.snooblish.daily.event.age <- table(floor(pop.reg[active.rows[which(is.snooblish & !is.snoob & who.died)], "age"]/365))
            M.nonsnoob.snooblish.events[names(nonsnoob.snooblish.daily.event.age)] <- M.nonsnoob.snooblish.events[names(nonsnoob.snooblish.daily.event.age)]+nonsnoob.snooblish.daily.event.age
        }
        
        
        
        
    }
    
    active.rows <- which(pop.reg[,"state"]==1)
    
    # emigration
    todays.n.emigrants <- daily.emigration(crude.emigration.rate, pop.data)
    names.who.left <- emigrant.picker(todays.n.emigrants, data=pop.reg)
    who.left <- rep(0, length(active.rows))
    names(who.left) <- active.rows
    who.left[as.character(names.who.left)] <- 1
    census.n.emigrants <- census.n.emigrants + sum(who.left)
    who.left <- as.logical(who.left)
    if(any(who.left)){
        census.emigration.list <- c(census.emigration.list, active.rows[who.left])
        pop.reg[active.rows[who.left], "state"] <- 2
        pop.reg[active.rows[who.left],"last.em"] <- pop.reg[active.rows[who.left], "age"]  # record their age at emigration
        if(event.logging==TRUE){
            left.list <- active.rows[who.left]
            for(i in 1:length(left.list)){
                event.log <- c(event.log, paste("Day ", day, ": ", ifelse(pop.reg[left.list[i], "male"]==1, "Male ", "Female "), left.list[i], " has left the population.", sep=""))
            }
        }
    }
    
    active.rows <- which(pop.reg[,"state"]==1)
    
    # individual change
    n.individual.changers <- rbinom(1, length(active.rows), individual.change.rate)
    
    if(n.individual.changers>0){
        who.changes <- sample(active.rows, n.individual.changers)
        active.row.t1s <- pop.reg[active.rows,"snoob"]
        for(i in 1:n.individual.changers){
            model.sample.t1s <- sample(active.row.t1s, social.sample.size, replace=T)
            sample.frq <- mean(model.sample.t1s)
            prob.converting.to.snoob <- (sample.frq^conformity.bias.snoob)/(sample.frq^conformity.bias.snoob + (1-sample.frq)^conformity.bias.snoob)
            pop.reg[who.changes[i],"snoob"] <- sample(c(1,0), 1, prob=c(prob.converting.to.snoob, 1-prob.converting.to.snoob))
        }
    }
    
    

    # bookkeeping to close out the day
    
    if(diagnostic==TRUE){
        readline("-- hit enter --")   # for diagnostics
        print(pop.reg[1:(min(which(is.na(pop.reg[,1])))-1),])  # for diagnostics
    }

    # everyone is a day older, including those who died?
    pop.reg[,"age"] <- as.numeric(pop.reg[,"age"]) + 1 
    
    # pregnant women are one day closer to birth
    pregnant <- active.rows[!is.na(pop.reg[active.rows,"counter"])]
    if(length(pregnant)>0){
        pop.reg[pregnant, "counter"] <- pop.reg[pregnant, "counter"] - 1
    }
    
    census.n.closing <- sum(pop.reg[,"state"]==1, na.rm=T)
    
    # record the census, if it is time
    if(day %% census.census.length == 0){
    
        census.rows <- c(active.rows, census.death.list, census.emigration.list)
        
        temp <- cbind(rep(day/census.census.length, length(census.rows)), census.rows, pop.reg[census.rows,])
        colnames(temp)[1] <- "census"
        colnames(temp)[2] <- "id"
        census.end.census <- rbind(census.end.census, temp)

        census.table <- as.matrix(rbind(census.table, c(census.n.opening, census.n.births, census.n.deaths, census.n.immigrants, census.n.emigrants, census.n.closing)))
        rownames(census.table) <- 1:dim(census.table)[1]
        
        census.n.opening <- 0
        census.n.closing <- 0
        census.n.births <- 0
        census.n.deaths <- 0
        census.n.immigrants <- 0
        census.n.emigrants <- 0
        census.death.list <- integer(0)
        census.emigration.list <- integer(0)
    }
    
    # check to see if the population register is running out of space, and increase accordingly
    if(day %% 30 == 0){
        if(sum(is.na(pop.reg[,"state"])) < 1000){
            last.entry <- dim(pop.reg)[1]
            pop.reg[(last.entry+1):(last.entry+1000),] <- rep(NA, dim(pop.reg)[2]*1000)
        }
    }
    
    # visualizing population growth
    if(day %% 365==0 & running.pop.graphs==TRUE){
        x.axis.day <- c(x.axis.day, day)
        n <- c(n, length(active.rows))
        n.snoob <- c(n.snoob, sum(pop.reg[active.rows,"snoob"]))
        n.eng <- c(n.eng, sum(pop.reg[active.rows, "language"]))
        par(mfrow=c(1,2))
        plot(x.axis.day, n, pch=20, xaxt="n", type="l", xlab="Years", ylim=c(0, max(n)))
        points(x.axis.day, n.snoob, pch=20, xaxt="n", type="l", xlab="Years", col="blue")
        points(x.axis.day, n.eng, pch=20, xaxt="n", type="l", xlab="Years", col="orange")
        abline(v=c(0, census.census.length*(1:floor(day/census.census.length))), lty=2, col="lightblue")
        stop.day <- max(x.axis.day) + (365-max(x.axis.day) %% 365)    
    #    quarterly.tks <- seq(0, stop.day, by=91.25)
    #    axis(1, at=quarterly.tks, labels=F)
        yearly.tks <- seq(0, stop.day, by=365)
        axis(1, at=yearly.tks, labels=0:(stop.day/365), lwd.ticks=1.5)
                
        age.pyramid(floor(pop.reg[active.rows,"age"]/365), pop.reg[active.rows, "male"], by=5)
    }
    
    day <- day + 1
    # if(day %% 365 == 0) print(day/365)
    
    should.simulate.day <- day < days.max & sum(is.na(pop.reg[,1])) != 0
    
    if(n.snoob[length(n.snoob)]/n[length(n)] > 0.99) should.simulate.day <- F
    
} # end day loop

####### Export Census Data #######

dir_init('./output')

my.data <- census.end.census
write.csv(my.data, "./output/raw_data.csv", row.names=FALSE)

print('simulation complete')

if(diagnostics){
    par(mfrow=c(2,1))
    event.vec <- M.snoob.nonsnooblish.events
    risk.vec <- M.snoob.nonsnooblish.at.risk
    expected.vec <- age.specific.mortality

    cat.names <- as.character(c(0, 1, 1, 1, 1, sort(rep(seq(5, 100, 5), 5))))

    barplot(event.vec[1:101]/risk.vec[1:101]*365*1000, space=0)
    points(0:104, expected.vec[cat.names], pch=20)
}

# Diagnostics for Mortality: 
# 1. basic mortality by age; how many events per possible events over the course fo the simulation...two vectors as long as each age last birthday...first vector adds every day the n.at.risk....second vector records n.events as they happen
# 2. snoob only: 
# M.snoob.at.risk
# M.snoob.events
# 3. nonsnoob only:
# M.nonsnoob.at.risk
# M.nonsnoob.events
# 4. snooblish only 
# M.snooblish.at.risk
# M.snooblish.events
# 5. nonsnooblish only
# M.nonsnooblish.at.risk
# M.nonsnooblish.events
# 6. snooblish only among snoobs
# M.snoob.snooblish.at.risk
# M.snoob.snooblish.events
# 7. nonsnooblish only among snoobs
# M.snoob.nonsnooblish.at.risk
# M.snoob.nonsnooblish.events
# 8. snooblish only among nonsnoobs
# M.nonsnoob.snooblish.at.risk
# M.nonsnoob.snooblish.events
# 9. nonsnooblish only among nonsnoobs
# M.nonsnoob.nonsnooblish.at.risk
# M.nonsnoob.nonsnooblish.events

# once the simulation is completed, these data vectors will tell me exactly what I want to know for mortality...its the only way to know for certain what's happening...
# if I had a "death matrix" I could compute the event counts manually....but if I'm keeping track of at risk then I might as well do it that way....blah

# fertility would be similar....
# the event here would be birthing

# Diagnostics for Births: 
# 1. basic mortality by age; how many events per possible events over the course fo the simulation...two vectors as long as each age last birthday...first vector adds every day the n.at.risk....second vector records n.events as they happen
# 2. snoob only: 
# B.snoob.at.risk
# B.snoob.events/365
# 3. nonsnoob only:
# B.nonsnoob.at.risk
# B.nonsnoob.events/365
# 4. snooblish only 
# B.snooblish.at.risk
# B.snooblish.events/365
# 5. nonsnooblish only
# B.nonsnooblish.at.risk
# B.nonsnooblish.events/365
# 6. snooblish only among snoobs
# B.snoob.snooblish.at.risk
# B.snoob.snooblish.events/365
# 7. nonsnooblish only among snoobs
# B.snoob.nonsnooblish.at.risk
# B.snoob.nonsnooblish.events/365
# 8. snooblish only among nonsnoobs
# B.nonsnoob.snooblish.at.risk
# B.nonsnoob.snooblish.events/365
# 9. nonsnooblish only among nonsnoobs
# B.nonsnoob.nonsnooblish.at.risk
# B.nonsnoob.nonsnooblish.events/365

# Diagnostics for Mating Assortment:  
# mated <- active.rows[!is.na(pop.reg[active.rows, "mate"])]
# cbind(pop.reg[mated, "snoob"], pop.reg[pop.reg[mated,"mate"], "snoob"])
# mean(as.numeric(pop.reg[mated, "snoob"]==pop.reg[pop.reg[mated,"mate"], "snoob"]))

# Diagnostics for Transmission Bias: 
# for each census, record midparent value for snoob and snooblish, and record kid value...then record the difference....

# Diagnostics for Individual Change:
# for each census, record the number of social learners and the average direction of the social learning.  Compare to the frq of snoobism and snooblish.  








# Exploring what the control panel parameters imply 

par(mfrow=c(1,3))

# social learning visualization
sample.frqs <- seq(0,1,0.001)
prob.converting.to.snoob <- (sample.frqs^conformity.bias.snoob)/(sample.frqs^conformity.bias.snoob + (1-sample.frqs)^conformity.bias.snoob)
plot(sample.frqs, prob.converting.to.snoob, type="l")
abline(0,1,lty=2)

# visualizing the mortality situation
x.points <- c(0, 1, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 100)
plot(x.points, age.specific.mortality, log="y", type="l", las=1, yaxt="n", ylim=c(0.01, 1000), ylab="Annual deaths per 1000 people", xlab="Age Last Birthday", main="Age-specific Mortality")
ticks <- 0:20*100
axis(1, at=ticks)
axis(2, at=c(0.01*(10^(0:5))), labels=c(0.01, 0.1, 1, 10, 100, 100), las=1)
baseline <- (age.specific.mortality/1000)/365
alpha = log(baseline/(1- baseline))

# try messing with these four parameters
male <- 0
snoob <- 1
daily.pr.death <- (1-(1/(1+exp(alpha + mortality.bias.male*male + mortality.bias.snoob*snoob))))
new.age.specific.mortality <- (daily.pr.death*1000)*365
points(x.points, new.age.specific.mortality, col="blue", type="l")

# visualizing fertility
baseline <- (age.specific.fertility/1000)/365
alpha = log(baseline/(1- baseline))
snoob <- 1
daily.pr.conception <- (1-(1/(1+exp(alpha + fertility.bias.snoob*snoob))))
new.age.specific.fertility <- (daily.pr.conception*1000)*365
x.points <- c(15, 20, 25, 30, 35, 40, 45)
plot(x.points, new.age.specific.fertility, type="l", col="blue", main=c("Age-specific fertility for Snoobs (blue) and non-Snoobs (black)"))
points(x.points, age.specific.fertility, col="black", type="l")


# TFR 
sum(age.specific.fertility*5)
# new TFR
sum(new.age.specific.fertility*5)


# so ends the pre-simulation EDA





# source of the immigrant effect: some of the initial population is not old enough to appear in the age pyramid...when they eventually DO appear, it is several censuses in and so they are counted as immigrants....this is probably contributing to the bugs...I need to ensure that the first 15 years are ignored, I guess.  Alternatively, I could just remove the kids under 15 from the initial population...don't kill them, though.








# For the simulator: 
# 1. get the simulation to record a list of all mating pairs to compare snoob status
# 2. ok, here's a question...if it is drawing from men at random, then how can we expect the age-specific fertility to apply to men then?  Does it at all, or is it only applying to women?
# run the ancestors function on the population register! 

# the code should not refer to the census as "census" - it's the census dammit!  And there shouldn't be a census 0...it should start at census 1...

# calculate age-specific mortality..that is, number of individuals who did die of each age category out of those who could have died and then vary by subpopulation.
# I need to check that mating pairs are as we'd expect given outgroup exclusion re snoob status.  THis seems like an in-simulation thing to check...  
# we need to establish what the mortality and fertility profiles are like for subpopulations.  The prepped data tells us if they had any kids intercensus, but we need to know how to credit that.  The censuses are spaced every five years, so they switched between age categories intercensus too.  I need to know mom and dad age at the exact moment the child was conceived to calculate this properly.  That would be something the simulation itself would have to have.  
# things the diagnostics have to do:
# is the transmission bias working properly?
# is the mate selection function working properly?
# etc. etc.

# in the simulator, when somenone dies, all their values that get stored in the census output should go to NA 

# ah, here's an interesting question: when i specify a variable as the grouping variable and it looks at someone who died, does it correctly look at their grouping variable from the last census they were alvie for or incorrectly look at the one at the moment they died?  That's a serios bug.  The simplest solution: use the id and state==0 to see who died, and from that moment on only look at the previous census entry for them.   What you should really do is NA all of their values when they are recorded as dead: age, age.cat, height, snoob, language, and risk.  

# some of these must be recorded as events when they happen...other things can be performed after the simulation ends, on the population register.  

# the million dollar question: 
# among nonsnoobs, is there in any evidence that snooblish speakers enjoy different mortality or fecundity profiles?  Likewise, among snoobs, is there evidence that non-snooblish speakers are in any way different from snooblish speakers?  

# within and among each age category....risk score?  

# within and between snoobs/nonsnoobs...


####
# the only way I can see how this could be done is by making a record when each event occurs...otherwise we won't know exactly who, when, and what they were and what the population looked like...
####























