
logistic <- function(x) (1 - (1/(1 + exp(x))))

pr_snoob <- function(p, beta) p^beta / (p^beta + (1 - p)^beta)

conception <- function(asf, reg, fecund.women.rows){
  
  age_bins <- reg$age_bin[fecund.women.rows]
  snoob <- reg$snoob[fecund.women.rows]
  
  baseline <- (asf[as.character(age_bins)]/1000)/365
  alpha = log(baseline/(1- baseline))
  
  daily.pr.conception <- logistic(alpha + fertility.bias.snoob * snoob)
  got.pregnant <- rbinom(length(fecund.women.rows), 1, daily.pr.conception)
  return(got.pregnant)

}

grim.reaper <- function(asm, reg, mortality.bias.male, mortality.bias.snoob) {

  alive <- which(reg$is_alive)

  baseline <- (asm[as.character(reg$age_bin[alive])]/1000)/365
  alpha = log(baseline/(1- baseline))

  daily.pr.death <- logistic(alpha + mortality.bias.male*reg$male[alive] + mortality.bias.snoob*reg$snoob[alive])
  died <- rbinom(length(alive), 1, daily.pr.death)
  
  who_died <- alive[which(died == 1)]
  return(who_died)
}

age.binner <- function(age_days){

  age.years <- floor(age_days/365)
  # 5-year age categories
  age_bins <- age.years - age.years %% 5
  age_bins[which(age_bins == 0)] <- 1
  age_bins[which(age_bins < 0)] <- 0 # ???
  age_bins[which(age_bins > 100)] <- 100

  return(age_bins)

}


baby.maker <- function(mom, dad, reg){

  out <- data.frame(
    mom = mom,
    dad = dad,
    is_alive = TRUE,
    age = 0,
    age_days = 0,
    age_bin = 0,
    male = rbinom(1, 1, prob = 0.5)
  )

  # snoob inheritance
  mom.snoob <- reg$snoob[mom]
  dad.snoob <- reg$snoob[dad]
  midparent.snoob <- mean(c(mom.snoob, dad.snoob))
  if (midparent.snoob == 1) midparent.snoob <- 0.9999
  if (midparent.snoob == 0) midparent.snoob <- 0.0001
  # this doesn't seem right btw
  baseline.trans.bias <- midparent.snoob
  alpha <- log(baseline.trans.bias / (1 - baseline.trans.bias))
  prob.kid.is.snoob <- logistic(alpha + transmission.bias.snoob)
  out$snoob <- rbinom(1, 1, prob = prob.kid.is.snoob)

  return(out)
}


mate.finder <- function(mom, reg) {
  # this seems the most inefficient of the simulator's code...it has to run mom by mom, which is quite a lot of them...
  
  current.mate <- NA
  
  # if the last mate exists, is alive and within the age range, he's your man
  if (!is.na(reg$mate[mom])) {
    if (reg$is_alive[reg$mate[mom]] & reg$age_days[reg$mate[mom]] > 15*365 & reg$age_days[reg$mate[mom]] < 60*365) {
      current.mate <- reg$mate[mom]
    }
  } else {
    available_men <- which(reg$male == 1 & reg$is_alive & reg$age_days >= 14 * 365 & reg$age_days < 60 * 365)
    if (length(available_men) > 0) {
      available_men.traits <- reg$snoob[available_men]
      n.men <- length(available_men)
      prob.choose.snoob <- mate.similarity.bias.snoob*reg$snoob[mom] + (1-reg$snoob[mom])*(1-mate.similarity.bias.snoob)
      prob.choose.nonsnoob <- 1 - prob.choose.snoob
      prob.is.chosen <- rep(0.5, n.men)
      prob.is.chosen[available_men.traits==1] <- prob.choose.snoob
      prob.is.chosen[available_men.traits==0] <- prob.choose.nonsnoob
      current.mate <- sample_safe(available_men, 1, prob=prob.is.chosen)
    }
  }
  current.mate
}
