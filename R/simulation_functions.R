
# workhorse functions

conception <- function(asf, reg, fecund.women.rows){
  
  age.cats <- reg$age.cat[fecund.women.rows]
  snoob <- reg$snoob[fecund.women.rows]
  
  baseline <- (asf[as.character(age.cats)]/1000)/365
  alpha = log(baseline/(1- baseline))
  
  daily.pr.conception <- (1 - (1 / (1 + exp(alpha + fertility.bias.snoob * snoob))))
  got.pregnant <- rbinom(length(fecund.women.rows), 1, daily.pr.conception)
  return(got.pregnant)

}

grim.reaper <- function(asm, reg, mortality.bias.male, mortality.bias.snoob){

  male <- reg[active.rows,"male"]
  snoob <- reg[active.rows,"snoob"]
  age.cats <- reg[active.rows, "age.cat"]
  
  baseline <- (asm[as.character(age.cats)]/1000)/365
  alpha = log(baseline/(1- baseline))

  daily.pr.death <- (1-(1/(1+exp(alpha + mortality.bias.male*male + mortality.bias.snoob*snoob))))
  died <- rbinom(length(active.rows), 1, daily.pr.death)
  return(died)
  
}

age.binner <- function(age.vec){

  zeros <- which(age.vec < 365)
  hundreds <- which(age.vec >= 36500)

  age.years <- floor(age.vec/365)
  age.cats <- age.years - age.years %% 5
  age.cats[age.cats==0] <- 1
  age.cats[zeros] <- 0
  age.cats[hundreds] <- 100

  return(age.cats)

}


baby.maker <- function(mom.id, dad.id, data){
  state <- 1
  age <- 0
  age.cat <- 0
  died <- NA
  male <- rbinom(1,1,prob=0.5)
  mate <- NA
  counter <- NA
  
  # snoob inheritance
  mom.snoob <- data[mom.id, "snoob"]
  dad.snoob <- data[dad.id, "snoob"]
  midparent.snoob <- mean(c(mom.snoob, dad.snoob))  # GOD DAMN that was a hard bug to solve...
  if(midparent.snoob == 1) midparent.snoob <- 0.9999
  if(midparent.snoob == 0) midparent.snoob <- 0.0001
  baseline.trans.bias <- midparent.snoob
  alpha <- log(baseline.trans.bias/(1-baseline.trans.bias))
  prob.kid.is.snoob <- exp(alpha + transmission.bias.snoob)/(1+ exp(alpha + transmission.bias.snoob))
  snoob <- rbinom(1,1,prob=prob.kid.is.snoob)

  newborn <- c(state, mom.id, dad.id, age, age.cat, died, male, mate, counter, snoob)
  return(newborn)
}


mate.finder <- function(mom.id, reg){  # this seems the most inefficient of the simulator's code...it has to run mom by mom, which is quite a lot of them...
  
  current.mate <- NA
  
  moms.snoob <- reg[mom.id, "snoob"]
  
  last.mate <- reg[mom.id, "mate"]
  if(!is.na(last.mate)){
    if(reg[last.mate, "state"] == 1 & reg[last.mate, "age"] > 15*365 & reg[last.mate, "age"] < 60*365){
      current.mate <- last.mate
    }
  }
  # if the last mate exists, is alive and within the age range, he's your man

  if(is.na(current.mate)){  
    # "available" means within age range, alive, male
    # male.mating.bias <- 1
    
    available.men <- active.rows[which( reg[active.rows, "male"]==1 & reg[active.rows, "age"] >= reproductive.min.age*365 & reg[active.rows, "age"] < 60*365 )]
          
    if(length(available.men)>0){
      current.mate <- available.men[1]
      if(length(available.men) > 1){
        available.men.traits <- reg[available.men,"snoob"]
        n.men <- length(available.men)
      
        prob.choose.snoob <- mate.similarity.bias.snoob*moms.snoob + (1-moms.snoob)*(1-mate.similarity.bias.snoob)
        prob.choose.nonsnoob <- 1 - prob.choose.snoob
  
        prob.is.chosen <- rep(0.5, n.men)
        prob.is.chosen[available.men.traits==1] <- prob.choose.snoob
        prob.is.chosen[available.men.traits==0] <- prob.choose.nonsnoob
                
        current.mate <- sample_safe(available.men, 1, prob=prob.is.chosen)
      }
    } else {
      current.mate <- NA
    }
  }
  current.mate  

}



# height function: 

height.fun <- function(age, h.gene){
  if(any(age>100)) age <- age/365
  heights <- h.gene*exp(-3 + .2*age)/(1 + exp(-3+.2*age))
  heights
}
