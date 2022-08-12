
rm(list = ls())

source("project_support.R")

cat("define simulation parameters\n")

# n_ind_init, n_sim_years, census_interval_years all defined in project_support.R
n_sim_days <- 365 * n_sim_years
census_interval <- 365 * census_interval_years

show_running_graphs <- FALSE
log_events <- FALSE

# Vital Population Rates
asf <- c(40.47, 102.50, 115.92, 96.10, 46.36, 9.12, 0.58)
# expected births per person within five-year age classes from 15 to 50
asm <- c(6.54, 0.29, 0.14, 0.18, 0.64, 0.91, 0.90, 1.06, 1.53, 2.31, 3.41, 4.93, 7.42, 11.50, 17.80, 27.71, 43.50, 69.58, 110.56, 174.77, 276.56, 438.92)
names(asf) <- c(15, 20, 25, 30, 35, 40, 45)
names(asm) <- c(0, 1, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 100)

# expected deaths per person within five-year age classes from 0 to 100+
daily_update_rate <- 0.00006  # fraction of individuals who update each day
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
plot(x.points, asm, log="y", type="l", las=1, yaxt="n", ylim=c(0.01, 1000), ylab="Annual deaths per 1000 people", xlab="Age Last Birthday", main="Age-specific Mortality")
ticks <- 0:20*100
axis(1, at=ticks)
axis(2, at=c(0.01*(10^(0:5))), labels=c(0.01, 0.1, 1, 10, 100, 100), las=1)
baseline <- (asm/1000)/365
alpha = log(baseline/(1- baseline))

# try messing with these four parameters
male <- 0
snoob <- 1
daily.pr.death <- (1-(1/(1+exp(alpha + mortality.bias.male*male + mortality.bias.snoob*snoob))))
new.asm <- (daily.pr.death*1000)*365
points(x.points, new.asm, col="blue", type="l")

# visualizing fertility
baseline <- (asf/1000)/365
alpha <- log(baseline/(1- baseline))
snoob <- 1
daily.pr.conception <- (1-(1/(1+exp(alpha + fertility.bias.snoob*snoob))))
new.asf <- (daily.pr.conception*1000)*365
x.points <- c(15, 20, 25, 30, 35, 40, 45)
plot(x.points, new.asf, type="l", col="blue", main=c("Age-specific fertility for Snoobs (blue) and non-Snoobs (black)"))
points(x.points, asf, col="black", type="l")

dev.off()



cat("initalize simulation\n")

# initial age structure, for calendar ages 0 to 99
age_weights <- c(57, 74, 68, 77, 65, 75, 79, 83, 66, 60, 80, 72, 59, 87, 78, 73, 58, 73, 84, 67, 83, 76, 71, 74, 72, 62, 91, 64, 63, 76, 80, 79, 71, 73, 76, 72, 70, 82, 80, 70, 82, 72, 78, 75, 75, 83, 73, 70, 79, 71, 77, 81, 80, 71, 64, 70, 62, 65, 73, 76, 70, 73, 84, 74 ,64 ,63 ,76 ,72 ,55 ,81 ,78 ,74 ,68 ,69 ,71 ,63 ,56, 50, 55, 12, 41, 38, 56, 30, 27, 22, 38, 30,  8, 22, 11 ,12, 12, 10,  4,  9,  8,  3,  2,  4)

reg <- data.frame(
  is_alive = rep(TRUE, n_ind_init),
  mom = NA,
  dad = NA,
  age = sample_safe(0:99, n_ind_init, replace = TRUE, prob = age_weights),
  age_at_death = NA,
  male = rbinom(n_ind_init, 1, 0.5),
  mate = NA,
  birth_counter = NA,
  snoob = rbinom(n_ind_init, 1, p=frq_snoob_init)
)

reg$age_days <- reg$age * 365 + sample_safe(0:364, n_ind_init, replace = TRUE)
reg$age_bin <- age.binner(reg$age_days)
reg$day_of_birth <- 0 - reg$age_days # simulation starts on day 0

stopifnot(all(c("is_alive", "mom", "dad", "age", "age_days", "age_bin", "age_at_death", "male",  "mate", "birth_counter", "snoob") %in% colnames(reg)))

if (log_events) event_log <- character(0)

# necessary for the figures
if (show_running_graphs) {
  x.axis.day <- 1
  n_ind <- n_ind_init
  n_snoob <- sum(reg$snoob)
}



cat("begin simulation daily loop\n")

cat("census", 1, "recorded on day", 0, "\n")
cens <- reg
cens$census <- 1
cens$day <- 0
cens$id <- seq_len(nrow(cens))

for (day in seq_len(n_sim_days)) {

  alive <- which(reg$is_alive)

  if (!any(reg$is_alive)) break()

  reg$age_bin[alive] <- age.binner(reg$age_days[alive])
  reg$age[alive] <- floor(reg$age_days[alive] / 365)

  # births
  birth_today <- which(reg$birth_counter == 0 & reg$is_alive)

  if (length(birth_today) > 0) {
    for (i in 1:length(birth_today)) {
      mom <- birth_today[i]
      dad <- reg$mate[birth_today[i]] # note that dads can die before kid is born!
      baby <- baby.maker(mom, dad, reg)
      baby$day_of_birth <- day
      reg <- bind_rows(reg, baby)
      if (log_events) {
        event_log <- c(event_log, paste("Day ", day, ": ", ifelse(baby$male == 1, "Male ", "Female "), nrow(reg), " has been born to mom ", baby$mom , " and dad ", baby$dad, ".", sep=""))
      }
      reg$birth_counter[mom] <- NA
    }
    alive <- which(reg$is_alive)
  }

  # mating
  # first, decide which women are going to conceive according to the age-specific fertility schedule...
  fecund <- which(reg$is_alive & reg$male == 0 & is.na(reg$birth_counter) & reg$age_bin >= 15 & reg$age_bin < 50)

  if (length(fecund) > 0) {
    will.conceive <- conception(asf, reg, fecund)
    # vector of length fecund of 0's and 1's
    female.maters <- fecund[which(will.conceive == 1)]

    # if any will conceive, those women will then attempt to select a partner...
    if (length(female.maters) > 0) {
      for (i in 1:length(female.maters)) {
        mom <- female.maters[i]
        has_mate <- !is.na(reg$mate[mom])
        dad <- mate.finder(mom, reg)
        reg$mate[mom] <- dad
        if (is.na(dad)) {
          if (log_events) {
            event_log <- c(event_log, paste("Day ", day, ": ", "Female ", mom, " tried to conceive but was unable to find a mate.", sep=""))
          }
        } else {
          reg$mate[dad] <- mom
          reg$birth_counter[mom] <- round(rnorm(1, 280, 5))
          if (log_events) {
            event_log <- c(event_log, paste("Day ", day, ": ", "Female ", mom, " has become pregnant from male ", dad, ".", sep=""))
          }
        }
      }
    }
  }

  # mortality

  died <- grim.reaper(asm, reg, mortality.bias.male, mortality.bias.snoob)
  # returns a vector of the rows (i.e. names) of the dead

  if (length(died) > 0) {
    reg$is_alive[died] <- FALSE
    reg$age_at_death[died] <- reg$age_days[died]  # record their age at death
    if (log_events) {
      for (i in 1:length(died)) {
        event_log <- c(event_log, paste("Day ", day, ": ", ifelse(reg$male[died[i]]==1, "Male ", "Female "), died[i], " has died.", sep=""))
      }
    }
    alive <- which(reg$is_alive)
  }
  
  # individual change
  n_updaters <- rbinom(1, length(alive), daily_update_rate)
  if (n_updaters > 0) {
    updaters <- sample_safe(alive, n_updaters)
    available_models <- reg$snoob[alive]
    for (i in 1:n_updaters) {
      my_models <- sample_safe(available_models, social.sample.size, replace = TRUE)
      sample.frq <- mean(my_models)
      # hmm, so people don't update based on their existing phenotypes!
      prob.update.to.snoob <- pr_snoob(sample.frq, conformity.bias.snoob)
      reg$snoob[updaters[i]] <- rbinom(1, 1, prob.update.to.snoob)
    }
  }

  # everyone alive is a day older
  reg$age_days[alive] <- reg$age_days[alive] + 1 
  
  # pregnant women are one day closer to birth
  pregnant <- which(reg$is_alive & !is.na(reg$birth_counter))
  if (length(pregnant) > 0) reg$birth_counter[pregnant] <- reg$birth_counter[pregnant] - 1
  # women who were at 0 before this section should have given birth and reset counter to NA
  stopifnot(all(is.na(reg$birth_counter) | reg$birth_counter >= 0))
  stopifnot(all(is.na(reg$birth_counter[which(reg$male == 1)])))

  # record the census
  if (day %% census_interval == 0) {
    census_number <- as.integer(day/census_interval) + 1
    cat("census", census_number, "recorded on day", day, "\n")
    add <- reg[alive,]
    add$census <- census_number
    add$id <- seq_len(nrow(reg))[alive]
    add$day <- day
    cens <- bind_rows(cens, add)
  }

  # visualizing population growth
  if (show_running_graphs) {
    if (day %% 365 == 0) {
      x.axis.day <- c(x.axis.day, day)
      n_ind <- c(n_ind, length(alive))
      n_snoob <- c(n_snoob, sum(reg$snoob[alive]))
      par(mfrow=c(1,2))
      plot(x.axis.day, n_ind, pch=20, xaxt="n", type="l", xlab="Years", ylim=c(0, max(n_ind)))
      points(x.axis.day, n_snoob, pch=20, xaxt="n", type="l", xlab="Years", col="blue")
      abline(v=c(0, census_interval*(1:floor(day/census_interval))), lty=2, col="lightblue")
      stop.day <- max(x.axis.day) + (365-max(x.axis.day) %% 365)
      yearly.tks <- seq(0, stop.day, by=365)
      axis(1, at=yearly.tks, labels=0:(stop.day/365), lwd.ticks=1.5)
      age.pyramid(floor(reg$age[alive]), reg$male[alive], by=5)
    }
  }

  if(day %% 365 == 0) cat("year", floor(day/365), "\n")

}

stopifnot(!any(diff(reg$day_of_birth[reg$day_of_birth > 0] < 0)))

reg$day_of_death <- reg$day_of_birth + reg$age_at_death

stopifnot(all(reg$day_of_birth < reg$day_of_death[reg$mom], na.rm = TRUE))

# dads can be dead when their kids are born - does that complicate the calculations?
diff_dad <- reg$day_of_birth - reg$day_of_death[reg$dad]
tar <- which(diff_dad > 0)
if (length(tar) > 0) stopifnot(any(diff_dad[tar] < (2 * 365)))

write.csv(cens, "censuses.csv", row.names=FALSE)
