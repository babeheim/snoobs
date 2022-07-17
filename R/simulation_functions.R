
# workhorse functions

daily.immigration <- function(crude.immigration.rate.f, data){
  daily.immigration.rate.per.1000 <- (crude.immigration.rate.f/1000)/365
  pop.size <- length(active.rows)
  mean.daily.n.immigrants <- daily.immigration.rate.per.1000*pop.size 
  todays.n.immigrants <- rpois(1, mean.daily.n.immigrants)    
  unoccupied.rows <- sum(is.na(data[,1]))
  return(todays.n.immigrants)
}

daily.emigration <- function(crude.emigration.rate.f, data){
  daily.emigration.rate.per.1000 <- (crude.emigration.rate.f/1000)/365
  pop.size <- length(active.rows)
  mean.daily.n.emigrants <- daily.emigration.rate.per.1000*pop.size
  todays.n.emigrants <- rpois(1, mean.daily.n.emigrants)
  return(todays.n.emigrants)    
}

emigrant.picker <- function(todays.n.emigrants.f, data){
  baseline.prob.emig <- 1/length(active.rows)
  t1s <- pop.reg[active.rows, "snoob"]
  alpha <- log(baseline.prob.emig/(1- baseline.prob.emig))
  emig.prob.vec <- exp(alpha + emigration.bias.snoob*t1s)/(1+exp(alpha + emigration.bias.snoob*t1s))
  emig.prob.vec <- emig.prob.vec/sum(emig.prob.vec)
  emigrants <- sample(active.rows, todays.n.emigrants.f, prob=emig.prob.vec)
  return(emigrants)
}

conception <- function(asf, data){
  
  age.cats <- data[fecund.women.rows, "age.cat"]
  snoob <- data[fecund.women.rows,"snoob"]
  
  baseline <- (asf[as.character(age.cats)]/1000)/365
  alpha = log(baseline/(1- baseline))
  
  daily.pr.conception <- (1-(1/(1+exp(alpha + fertility.bias.snoob*snoob))))
  got.pregnant <- rbinom(length(fecund.women.rows), 1, daily.pr.conception)
  return(got.pregnant)

}

grim.reaper <- function(asm, data){

  male <- data[active.rows,"male"]
  snoob <- data[active.rows,"snoob"]
  age.cats <- data[active.rows, "age.cat"]
  
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



immigrant.maker <- function(data){ 
  state <- 1
  mom <- 0
  dad <- 0
  age <- sample(1:(365*50),1)
  age.cat <- age.binner(age)
  last.im <- age
  last.em <- NA
  died <- NA
  male <- rbinom(1,1,prob=0.5)
  mate <- NA
  counter <- NA
  
  h.gene <- sample( 80:120 , 1)
  height <- height.fun(age=age, h.gene=h.gene) 
  
  baseline.immigration <- mean(data[active.rows, "snoob"])
  alpha <- log(baseline.immigration/(1-baseline.immigration))
  
  snoob <- rbinom(1, 1,prob=(exp(alpha + immigration.bias.snoob)/(1+exp(alpha + immigration.bias.snoob))))
  
  language <- rbinom(1, 1, mean(data[active.rows, "language"]))
  
  aiy <- age/365
  a <- 0.18
  if(snoob==1) a <- a*snoob.less.risky.factor
  b <- 0.07
  risk <- a*aiy/exp(b*aiy)
  
  normal.traits <- c(snoob, language, risk)
  
  trait.vec <- c(h.gene, height, normal.traits)
    
  immigrant <- c(state, mom, dad, age, age.cat, last.im, last.em, died, male, mate, counter, trait.vec)
  immigrant
}


risk.setter <- function(data){
  aiy <- pop.reg[active.rows,"age"]/365
  snoob <- pop.reg[active.rows,"snoob"]
  a <- 0.18
  b <- 0.07
  a.vec <- rep(a, length(active.rows))
  a.vec[snoob==1] <- snoob.less.risky.factor*a
  risk <- a.vec*aiy/exp(b*aiy)
  return(risk)
}

# I should set the risk scoring to be (a) stochastic and (b) peaking at the same time as RS, which is around 25 yrs old.  

# a <- 0.18
# b <- 0.07
# curve(100*a*x/exp(b*x), from=0, to=90, add=T)
# abline(v=1/b, lty=2, col="gray")

baby.maker <- function(mom.id, dad.id, data){
  state <- 1
  age <- 0
  age.cat <- 0
  last.im <- NA
  last.em <- NA
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
  
  # english inheritance
  mom.english <- data[mom.id, "language"]
  dad.english <- data[dad.id, "language"]
  parent.frq.english <- mean(c(mom.english, dad.english))
  parent.alpha <- 0.99
  pop.frq.english <- mean(pop.reg[active.rows,"language"])
  pr.acquire.english <- parent.alpha*parent.frq.english + (1-parent.alpha)*pop.frq.english
  language <- rbinom(1, 1, pr.acquire.english)
  
  risk <- 0
  
  normal.traits <- c(snoob, language, risk)
  
  h.gene <- sample( c(data[mom.id, "h.gene"], data[mom.id, "h.gene"]) , 1)
  height <- height.fun(age=age, h.gene=h.gene) 
  
  trait.vec <- c(h.gene, height, normal.traits)
  
  
  newborn <- c(state, mom.id, dad.id, age, age.cat, last.im, last.em, died, male, mate, counter, trait.vec)
  newborn
}


mate.finder <- function(mom.id, data){  # this seems the most inefficient of the simulator's code...it has to run mom by mom, which is quite a lot of them...
  
  current.mate <- NA
  
  moms.snoob <- data[mom.id, "snoob"]
  
  last.mate <- data[mom.id, "mate"]
  if(!is.na(last.mate)){
    if(data[last.mate, "state"] == 1 & data[last.mate, "age"] > 15*365 & data[last.mate, "age"] < 60*365){
      current.mate <- last.mate
    }
  }
  # if the last mate exists, is alive and within the age range, he's your man
  
  if(is.na(current.mate)){  
    # "available" means within age range, alive, male
    # male.mating.bias <- 1
    
    available.men <- active.rows[which( data[active.rows, "male"]==1 & data[active.rows, "age"] >= reproductive.min.age*365 & data[active.rows, "age"] < 60*365 )]
          
    if(length(available.men)>0){
      current.mate <- available.men[1]
      if(length(available.men) > 1){
        available.men.traits <- data[available.men,"snoob"]
        n.men <- length(available.men)
      
        prob.choose.snoob <- mate.similarity.bias.snoob*moms.snoob + (1-moms.snoob)*(1-mate.similarity.bias.snoob)
        prob.choose.nonsnoob <- 1 - prob.choose.snoob
  
        prob.is.chosen <- rep(0.5, n.men)
        prob.is.chosen[available.men.traits==1] <- prob.choose.snoob
        prob.is.chosen[available.men.traits==0] <- prob.choose.nonsnoob
                
        current.mate <- sample(available.men, 1, prob=prob.is.chosen)
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
