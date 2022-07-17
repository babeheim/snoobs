
cat("define simulation parameters\n")

# n_ind_init, n_sim_years, census_interval_years all defined in project_support.R
final_day <- 365 * n_sim_years
census_interval <- 365 * census_interval_years

show_running_graphs <- FALSE
log_events <- FALSE

# Vital Population Rates
reproductive.min.age <- 14 # in years
age.specific.fertility <- c(40.47, 102.50, 115.92, 96.10, 46.36, 9.12, 0.58)
# expected births per person within five-year age classes from 15 to 50
age.specific.mortality <- c(6.54, 0.29, 0.14, 0.18, 0.64, 0.91, 0.90, 1.06, 1.53, 2.31, 3.41, 4.93, 7.42, 11.50, 17.80, 27.71, 43.50, 69.58, 110.56, 174.77, 276.56, 438.92)
# expected deaths per person within five-year age classes from 0 to 100+
update_rate <- 0.00006  # fraction of individuals who update each day
social.sample.size <- 100  # number of individuals sampled when updating
mortality.bias.male <- 0.3 # logistic coefficient: positive values mean men die more often

# Mechanisms of Evoluton
mate.similarity.bias.snoob <- 0.9 # from 0 to 1; if 0.5, no assortment by snoob status.  If 1, perfect assortment along snoob status.  If 0, perfect reverse assortment along snoob status.
fertility.bias.snoob <- 0.8 # logistic coefficient; positive values mean snoobs have higher RS than nonsnoobs, all else being equal
mortality.bias.snoob <- -0.5 # logistic coefficient: positive values mean snoobs die more often.
transmission.bias.snoob <- 0.0   # logistic coefficient: if 0, kid will take snoob=1 with probability equal to midparent value.  If >0, kid will be more likely to take snoob=1 than would be expected from the midparent value.  If <0, kid will be LESS likely to take snoob=1 than would be expected from midparent values....to see real effects, this term has to be LARGE - magnitude 10 or so.
conformity.bias.snoob <- 1.3 # exponential coefficient: unbiased transmission if equals 1....values > 1 mean biased copying towards majority's trait, values < 1 mean biased copying against majority's trait
snoob.less.risky.factor <- 0.9 # varies between 0 and 1; 1 indicates snoobs are just as risky as non-snoobs, while 0 indicates they aren't risky at all

stopifnot(mate.similarity.bias.snoob >= 0 & mate.similarity.bias.snoob <= 1)


cat("create figure simulation_schedules\n")

png("figures/simulation_schedules.png", res = 300, height = 5, width = 15, units = "in")

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
alpha <- log(baseline/(1- baseline))
snoob <- 1
daily.pr.conception <- (1-(1/(1+exp(alpha + fertility.bias.snoob*snoob))))
new.age.specific.fertility <- (daily.pr.conception*1000)*365
x.points <- c(15, 20, 25, 30, 35, 40, 45)
plot(x.points, new.age.specific.fertility, type="l", col="blue", main=c("Age-specific fertility for Snoobs (blue) and non-Snoobs (black)"))
points(x.points, age.specific.fertility, col="black", type="l")

dev.off()

cat("initalize simulation\n")

# SnoobSim State Codes:
# 0 - Dead 
# 1 - Active 

reg.max <- n_ind_init*2  # this isn't really important anymore...the matrix will grow

day <- 1
x.axis.day <- 1

reg_traits <- c("state", "mom", "dad", "age", "age.cat", "died", "male",  "mate", "counter", "snoob")

state <- c(rep(1, n_ind_init), rep(NA, reg.max-n_ind_init))
mom <- c(rep(0, n_ind_init), rep(NA, reg.max-n_ind_init))

reg <- matrix(NA,ncol=length(reg_traits), nrow=reg.max)
reg <- as.data.frame(reg)
colnames(reg) <- reg_traits

reg$state[1:n_ind_init] <- 1
reg$mom[1:n_ind_init] <- 0
reg$dad[1:n_ind_init] <- 0

# initial age structure
ages <- 0:99
age_weights <- c(57, 74, 68, 77, 65, 75, 79, 83, 66, 60, 80, 72, 59, 87, 78, 73, 58, 73, 84, 67, 83, 76, 71, 74, 72, 62, 91, 64, 63, 76, 80, 79, 71, 73, 76, 72, 70, 82, 80, 70, 82, 72, 78, 75, 75, 83, 73, 70, 79, 71, 77, 81, 80, 71, 64, 70, 62, 65, 73, 76, 70, 73, 84, 74 ,64 ,63 ,76 ,72 ,55 ,81 ,78 ,74 ,68 ,69 ,71 ,63 ,56, 50, 55, 12, 41, 38, 56, 30, 27, 22, 38, 30,  8, 22, 11 ,12, 12, 10,  4,  9,  8,  3,  2,  4)

reg$age[1:n_ind_init] <- 365 * sample_safe(ages, n_ind_init, replace = TRUE, prob = age_weights) + sample_safe(0:364, n_ind_init, replace = TRUE)
reg$age.cat[1:n_ind_init] <- age.binner(reg$age[1:n_ind_init])
reg$male[1:n_ind_init] <- rep(c(1,0), n_ind_init/2)
reg$snoob[1:n_ind_init] <- rbinom(n_ind_init, 1, p=frq_snoob_init)

census.death.list <- integer(0)
if (log_events) event.log <- character(0)

disable.births <- FALSE # for diagnostic purposes
disable.deaths <- FALSE  # for diagnostic purposes
if (disable.births) age.specific.fertility <- rep(0, 7)
if (disable.deaths) age.specific.mortality <- rep(0, 22)
active.rows <- which(reg$state==1)

census_registers <- cbind(rep(0, length(active.rows)), active.rows, reg[active.rows,])
colnames(census_registers)[1] <- "census"
colnames(census_registers)[2] <- "id"

names(age.specific.fertility) <- c(15, 20, 25, 30, 35, 40, 45)
names(age.specific.mortality) <- c(0, 1, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 100)

n_ind <- length(active.rows)
n_snoob <- sum(reg$snoob[active.rows])

should.simulate.day <- day < final_day & sum(is.na(reg[,1])) != 0

cat("begin simulation daily loop\n")

while (should.simulate.day) {

  active.rows <- which(reg$state==1)

  reg$age.cat[active.rows] <- age.binner(reg$age[active.rows])

  # births
  giving.birth <- active.rows[which(!is.na(reg$counter[active.rows]) & reg$counter[active.rows]==0)]
  if(length(giving.birth)>0){
    for(i in 1:length(giving.birth)){
      mom <- giving.birth[i]
      dad <- reg$mate[giving.birth[i]]
      baby.surname <- min(which(is.na(reg[,1])))
      reg[baby.surname, ] <- baby.maker(mom, dad, reg)
      if (log_events) {
        event.log <- c(event.log, paste("Day ", day, ": ", ifelse(reg$male[baby.surname]==1, "Male ", "Female "), baby.surname, " has been born to female ", reg$mom[baby.surname] , " and male ", reg$dad[baby.surname], ".", sep=""))
      }
      reg$counter[mom] <- NA
    }
  }

  active.rows <- which(reg$state==1)

  # mating
  male <- reg[active.rows, "male"]
  not.pregnant <- is.na(reg[active.rows, "counter"])
  active.age.cats <- reg[active.rows, "age.cat"]
  fecund.women.rows <- active.rows[which(male == 0 & not.pregnant & active.age.cats >= 15 & active.age.cats < 50)]
  will.conceive <- conception(age.specific.fertility, reg, fecund.women.rows)
  female.maters <- fecund.women.rows[as.logical(will.conceive)]

  if (length(female.maters) > 0){
    for (i in 1:length(female.maters)) {
      mom <- female.maters[i]
      dad <- mate.finder(mom, reg)
      reg[mom,"mate"] <- dad
      if(is.na(dad)){
        if (log_events) {
          event.log <- c(event.log, paste("Day ", day, ": ", "Female ", mom, " tried to conceive but was unable to find a mate.", sep=""))
        }
      } else {
        reg[dad,"mate"] <- mom
        reg[mom, "counter"] <- round(rnorm(1,280,5))
        if (log_events) {
          event.log <- c(event.log, paste("Day ", day, ": ", "Female ", mom, " has become pregnant from male ", dad, ".", sep=""))
        }
      }
    }
  }

  # mortality
  who.died <- grim.reaper(age.specific.mortality, reg, mortality.bias.male, mortality.bias.snoob)
  who.died <- as.logical(who.died)

  if(any(who.died)){
    census.death.list <- c(census.death.list, active.rows[who.died])
    reg[active.rows[who.died], "state"] <- 0
    reg[active.rows[who.died],"died"] <- reg[active.rows[who.died], "age"]  # record their age at death
    if (log_events) {
      death.list <- active.rows[who.died]
      for(i in 1:length(death.list)){
        event.log <- c(event.log, paste("Day ", day, ": ", ifelse(reg[death.list[i], "male"]==1, "Male ", "Female "), death.list[i], " has died.", sep=""))
      }
    }
  }
  
  active.rows <- which(reg$state == 1)
  
  # individual change
  n_updaters <- rbinom(1, length(active.rows), update_rate)
  if (n_updaters > 0) {
    who.changes <- sample_safe(active.rows, n_updaters)
    active.row.t1s <- reg[active.rows,"snoob"]
    for(i in 1:n_updaters){
      model.sample.t1s <- sample_safe(active.row.t1s, social.sample.size, replace=T)
      sample.frq <- mean(model.sample.t1s)
      prob.converting.to.snoob <- (sample.frq^conformity.bias.snoob)/(sample.frq^conformity.bias.snoob + (1-sample.frq)^conformity.bias.snoob)
      reg[who.changes[i],"snoob"] <- sample_safe(c(1,0), 1, prob=c(prob.converting.to.snoob, 1-prob.converting.to.snoob))
    }
  }

  # bookkeeping to close out the day

  # everyone is a day older, including those who died?
  reg$age <- reg$age + 1 
  
  # pregnant women are one day closer to birth
  pregnant <- active.rows[!is.na(reg[active.rows,"counter"])]
  if(length(pregnant)>0){
    reg[pregnant, "counter"] <- reg[pregnant, "counter"] - 1
  }
  
  # record the census, if it is time
  if (day %% census_interval == 0) {
    census.rows <- c(active.rows, census.death.list)
    temp <- cbind(rep(day/census_interval, length(census.rows)), census.rows, reg[census.rows,])
    colnames(temp)[1] <- "census"
    colnames(temp)[2] <- "id"
    census_registers <- rbind(census_registers, temp)
    census.death.list <- integer(0)
  }
  
  # check to see if the population register is running out of space, and increase accordingly
  if (day %% 30 == 0){
    if (sum(is.na(reg$state)) < 1000) {
      last.entry <- dim(reg)[1]
      reg[(last.entry+1):(last.entry+1000),] <- rep(NA, dim(reg)[2]*1000)
    }
  }

  # visualizing population growth
  if (show_running_graphs) {
    if (day %% 365 == 0) {
      x.axis.day <- c(x.axis.day, day)
      n_ind <- c(n_ind, length(active.rows))
      n_snoob <- c(n_snoob, sum(reg[active.rows,"snoob"]))
      par(mfrow=c(1,2))
      plot(x.axis.day, n_ind, pch=20, xaxt="n", type="l", xlab="Years", ylim=c(0, max(n_ind)))
      points(x.axis.day, n_snoob, pch=20, xaxt="n", type="l", xlab="Years", col="blue")
      abline(v=c(0, census_interval*(1:floor(day/census_interval))), lty=2, col="lightblue")
      stop.day <- max(x.axis.day) + (365-max(x.axis.day) %% 365)
      yearly.tks <- seq(0, stop.day, by=365)
      axis(1, at=yearly.tks, labels=0:(stop.day/365), lwd.ticks=1.5)
      age.pyramid(floor(reg[active.rows,"age"]/365), reg[active.rows, "male"], by=5)
    }
  }

  day <- day + 1
  if(day %% 365 == 0) print(day/365)

  should.simulate.day <- (day < final_day) & (sum(is.na(reg[,1])) != 0)

}

write.csv(census_registers, "raw_data.csv", row.names=FALSE)
