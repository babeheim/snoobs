
if (scaffold) {
  rm(list = ls())
  source("../project_support.r")
}

dir_init("./temp")

# why does the prepped data have fewer rows than the raw data? because we have to drop babies who died before their first census - they were recorded in the raw data because that just records all intercensus deaths, but technically they don't "exist" in the analysis because they were never alive at a census, so they have to go.

# an interesting question: when you have a fit probability model, how do you plot expected outcomes?  You'd have to know the EV for the distribution, right?  I know a X ~ binomial's EV is p, and a X ~ normal's EV is mu, which allows me to plot linear Gaussian models and logistic models.  A poisson's EV is lambda.  Once I have a fit model, I can plot the EV against variation in some predictor variable like age.  But what if I don't know the model's EV??  I have a zero-inflated binomial and zero-inflated poisson in the mix...

raw.data <- read.csv("./inputs/raw_data.csv") # load the raw data
prepped.data <- read.csv("./inputs/prepped_data.csv") # also load in the prepped data (need both)

# I want to know how many generations are rerpseented in the data...
id.data <- raw.data[,"id"]
dad.data <- raw.data[,"dad"]
mom.data <- raw.data[,"mom"]
drop <- which(duplicated(id.data))

id.data <- id.data[-drop]
dad.data <- dad.data[-drop]
mom.data <- mom.data[-drop]

drop <- 1:1000

id.data <- id.data[-drop]
dad.data <- dad.data[-drop]
mom.data <- mom.data[-drop]




# make some custom functions


su <- function(data) (data - mean(data, na.rm=T))/sdd(data)

cen <- function(data){
	data - mean(data)
}

my.dzbinom <- function(x, size, prob, a, log=FALSE){
	output <- ifelse(x==0, a, 0) + (1-a)*choose(size,x)*prob^x*(1-prob)^(size-x)
	if(log==T) output <- log(output)
	output
}

my.dzpois <- function(x, lambda, a, log=FALSE){
	output <- ifelse(x==0, a, 0) + (1-a)*(exp(-lambda)*lambda^x)/factorial(x)
	if(log==T) output <- log(output)
	output
}

sdd <- function(data){
    x <- na.omit(data)
    sqrt(sum((x-mean(x))^2)/length(x))
}












# Mortality Analysis


# there's two ways to do it: predicted lifetime based on ...er...traits
# or more simply, prediction about whether you die or not next intercensus based on your traits NOW.  that's the one to do...

drop <- which(prepped.data[,"state"]==3)
living.data <- prepped.data[-drop,]

census.data <- living.data[,"census"]
census.list <- unique(living.data[,"census"])

# we need to know someone's status NEXT census
# we'll have to have data for alive at time t to prredict their death at t +1...so needs to be living data! 

snoob <- living.data[,"phi"]
will.die <- rep(0, nrow(living.data))
age <- integer(0)
male <- integer(0)

census.n <- tapply(snoob, census.data, length)
current.n <- census.n[as.character(census.data)]
census.phibar <- tapply(snoob, census.data, mean)
current.phibar <- census.phibar[as.character(census.data)]




for(i in 1:(length(census.list)-1)){
	census.names <- living.data[living.data[,"census"]==census.list[i], "id"]
	prepped.subset <- prepped.data[prepped.data[,"census"]==census.list[i+1],]
	will.die[which(census.data==census.list[i])] <- as.numeric(prepped.subset[match(census.names, prepped.subset[,"id"]),"state"]==3)
	raw.data.subset <- raw.data[raw.data[,"census"]==census.list[i],]
	male <- c(male, raw.data.subset[match(census.names, raw.data.subset[,"id"]),"male"])
	age <- c(age, raw.data.subset[match(census.names, raw.data.subset[,"id"]),"age"])
}

last.census.drop <- which(census.data==53)
snoob <- snoob[-last.census.drop]
will.die <- will.die[-last.census.drop]
age <- age[-last.census.drop]
male <- male[-last.census.drop]
current.n <- current.n[-last.census.drop]
current.phibar <- current.phibar[-last.census.drop]
census.data <- census.data[-last.census.drop]


age.yr <- age/365


m0 <- mle2(y ~ dbinom(1, prob=1/(1+exp(a))), data=list(y=will.die), start=list(a=0.1), trace=FALSE)

m1 <- mle2(y ~ dbinom(1, prob=1/(1+exp(p1 + p2*cen(age.yr)))) , data=list(y=will.die), start=list(p1=4.8, p2=0), trace=FALSE)

m2 <- mle2(y ~ dbinom(1, prob=1/(1+exp(p1 + p2*cen(age.yr) + p3*cen(age.yr)^2))) , data=list(y=will.die), start=list(p1=4.8, p2=0, p3=0), trace=FALSE)

m3 <- mle2(y ~ dbinom(1, prob=1/(1+exp(p1 + p2*cen(age.yr) + p3*cen(age.yr)^2 + p4*male))) , data=list(y=will.die), start=list(p1=4.8, p2=0, p3=0, p4=0), trace=FALSE)

m4 <- mle2(y ~ dbinom(1, prob=1/(1+exp(p1 + p2*cen(age.yr) + p3*cen(age.yr)^2 + p4*male + p5*snoob))) , data=list(y=will.die), start=list(p1=4.8, p2=0, p3=0, p4=0, p5=0), trace=FALSE)

m5 <- mle2(y ~ dbinom(1, prob=1/(1+exp(p1 + p2*cen(age.yr) + p3*cen(age.yr)^2 + p4*male + p5*snoob + p6*current.phibar))) , data=list(y=will.die), start=list(p1=4.8, p2=0, p3=0, p4=0, p5=0, p6=0), trace=FALSE)

m6 <- mle2(y ~ dbinom(1, prob=1/(1+exp(p1 + p2*cen(age.yr) + p3*cen(age.yr)^2 + p4*male + p5*snoob + p6*su(age.yr)^3))) , data=list(y=will.die), start=list(p1=4.8, p2=0, p3=0, p4=0, p5=0, p6=0), trace=FALSE)





m0z <- mle2(y ~ my.dzbinom(1, prob=1/(1+exp(p1)), a=alpha), data=list(y=will.die), start=list(p1=0.1, alpha=0.1), trace=FALSE)

m1z <- mle2(y ~ my.dzbinom(1, prob=1/(1+exp(p1 + p2*cen(age))), a=alpha) , data=list(y=will.die), start=list(p1=3.14, p2=0, alpha=0.1), trace=FALSE)

# m2z <- mle2(y ~ my.dzbinom(1, prob=1/(1+exp(p1 + p2*cen(age.yr) + p3*cen(age.yr)^2)), a=alpha) , data=list(y=will.die), start=list(p1=3.66, p2=0, p3=0, alpha=0.1), trace=FALSE)

m3z <- mle2(y ~ my.dzbinom(1, prob=1/(1+exp(p1 + p2*cen(age.yr) + p3*cen(age.yr)^2 + p4*male)), a=alpha) , data=list(y=will.die), start=list(p1=5.39, p2=0, p3=0, p4=0, alpha=0.1), trace=FALSE)

m4z <- mle2(y ~ my.dzbinom(1, prob=1/(1+exp(p1 + p2*cen(age.yr) + p3*cen(age.yr)^2 + p4*male + p5*snoob)), a=1/(1+exp(alpha))) , data=list(y=will.die), start=list(p1=4.8, p2=0, p3=0, p4=0, p5=0, alpha=0.1), trace=FALSE)

m5z <- mle2(y ~ my.dzbinom(1, prob=1/(1+exp(p1 + p2*cen(age.yr) + p3*cen(age.yr)^2 + p4*male + p5*snoob + p6*current.phibar)), a=1/(1+exp(alpha))) , data=list(y=will.die), start=list(p1=4.8, p2=0, p3=0, p4=0, p5=0, p6=0, alpha=0.1), trace=FALSE)

# AICtab(m0, m1, m2, m3, m4, m5, m0z, m1z, m2z, m3z, m4z, m5z, m6, weights=T)
# compare(m0, m1, m2, m3, m4, m5, m0z, m1z, m2z, m3z, m4z, m5z, m6, nobs=length(will.die))



age.yr.c <- floor(age.yr)
age.yr.c[age.yr.c > 100] <- 100



png("./temp/mortality.png", height=15, width=20, 
        units='cm', res=300)

par(mfrow=c(2,2))

s.pr.die.nums <- tapply(will.die[snoob==1], age.yr.c[snoob==1], mean)
s.pr.die <- rep(0, length(0:100))
names(s.pr.die) <- 0:100
s.pr.die[names(s.pr.die.nums)] <- s.pr.die.nums

ns.pr.die.nums <- tapply(will.die[snoob==0], age.yr.c[snoob==0], mean)
ns.pr.die <- rep(0, length(0:100))
names(ns.pr.die) <- 0:100
ns.pr.die[names(ns.pr.die.nums)] <- ns.pr.die.nums


p1 <- 4.673
p2 <- -0.036
p3 <- -0.001
p4 <- -0.088
p5 <- 0.559

plot(1000*ns.pr.die ~ names(ns.pr.die), las=1)
points(1000*s.pr.die ~ names(ns.pr.die), pch=20)
curve(1000*dbinom(1, 1, prob=1/(1+exp(p1 + p2*(x-33.81) + p3*(x-33.81)^2 + p4*0 + p5*0))), type="l", lty=2, from=0, to=100, add=T, col="maroon1")
curve(1000*dbinom(1, 1, prob=1/(1+exp(p1 + p2*(x-33.81) + p3*(x-33.81)^2 + p4*0 + p5*1))), type="l", from=0, to=100, add=T, col="maroon1")
curve(1000*dbinom(1, 1, prob=1/(1+exp(p1 + p2*(x-33.81) + p3*(x-33.81)^2 + p4*1 + p5*0))), type="l", lty=2, from=0, to=100, add=T, col="dodgerblue")
curve(1000*dbinom(1, 1, prob=1/(1+exp(p1 + p2*(x-33.81) + p3*(x-33.81)^2 + p4*1 + p5*1))), type="l", from=0, to=100, add=T, col="dodgerblue")

s.pr.die.nums <- tapply(will.die[snoob==1], age.yr.c[snoob==1], mean)
s.pr.die <- rep(0, length(0:100))
names(s.pr.die) <- 0:100
s.pr.die[names(s.pr.die.nums)] <- s.pr.die.nums

ns.pr.die.nums <- tapply(will.die[snoob==0], age.yr.c[snoob==0], mean)
ns.pr.die <- rep(0, length(0:100))
names(ns.pr.die) <- 0:100
ns.pr.die[names(ns.pr.die.nums)] <- ns.pr.die.nums


p1 <- 4.673
p2 <- -0.036
p3 <- -0.001
p4 <- -0.088
p5 <- 0.559

plot(1000*ns.pr.die ~ names(ns.pr.die), log="y", las=1, col="gray", ylim=c(.1, 1000))
points(1000*s.pr.die ~ names(ns.pr.die), pch=1)

x <- c(0:100)
points(1000*dbinom(1, 1, prob=1/(1+exp(p1 + p2*(x-33.81) + p3*(x-33.81)^2 + p4*0 + p5*0))), type="l", col="maroon1")
points(1000*dbinom(1, 1, prob=1/(1+exp(p1 + p2*(x-33.81) + p3*(x-33.81)^2 + p4*0 + p5*1))), type="l", col="maroon1")
points(1000*dbinom(1, 1, prob=1/(1+exp(p1 + p2*(x-33.81) + p3*(x-33.81)^2 + p4*1 + p5*0))), type="l", col="dodgerblue")
points(1000*dbinom(1, 1, prob=1/(1+exp(p1 + p2*(x-33.81) + p3*(x-33.81)^2 + p4*1 + p5*1))), type="l", col="dodgerblue")

human.ages <- 0:100
plot(cumprod(1-ns.pr.die))
points(cumprod(1-s.pr.die), pch=20)
fitted.s <- dbinom(1, 1, prob=1/(1+exp(p1 + p2*(human.ages-33.81) + p3*(human.ages-33.81)^2 + p4*0 + p5*1)))
fitted.ns <- dbinom(1, 1, prob=1/(1+exp(p1 + p2*(human.ages-33.81) + p3*(human.ages-33.81)^2 + p4*0 + p5*0)))
points(cumprod(1-fitted.s), type="l", col="maroon1")
points(cumprod(1-fitted.ns), type="l", col="maroon1")
fitted.s <- dbinom(1, 1, prob=1/(1+exp(p1 + p2*(human.ages-33.81) + p3*(human.ages-33.81)^2 + p4*1 + p5*1)))
fitted.ns <- dbinom(1, 1, prob=1/(1+exp(p1 + p2*(human.ages-33.81) + p3*(human.ages-33.81)^2 + p4*1 + p5*0)))
points(cumprod(1-fitted.s), type="l", col="dodgerblue")
points(cumprod(1-fitted.ns), type="l", col="dodgerblue")

# Snoobs tend to be younger! 

ns.age <- tapply(age[snoob==0], census.data[snoob==0], mean)
s.age <- tapply(age[snoob==1], census.data[snoob==1], mean)
plot(ns.age/365, type="l", ylim=c(25, 45), ylab="Age in Years", xlab="Census")
points(s.age/365, col="red", type="l")

dev.off()

















# Fertility Analysis

# Do I want to calculate "children even had" or just the kids had last intercensus?




# # basic models
# lambda=b0 + b1*age
# lambda=b0 + b1*age + b2*age^2
# lambda=b0 + b1*age + b2*age^2 + b3*snoob
# lambda=b0 + b1*age + b2*age^2 + b3*snoob + b4*male

# # with density-dependence
# lambda=(1-N/K)*(b0 + b1*age + b2*age^2 + b3*snoob + b4*male) # where N is current pop size
# lambda=exp(-alpha*N)*(b0 + b1*age + b2*age^2 + b3*snoob + b4*male) # where N is current pop size

# # with frequency-dependence
# lambda=b0 + b1*age + b5*phibar
# lambda=b0 + b1*age + b2*age^2 + b5*phibar
# lambda=b0 + b1*age + b2*age^2 + b3*snoob + b5*phibar
# lambda=b0 + b1*age + b2*age^2 + b3*snoob + b4*male + b5*phibar





# Kids will have next intercensus:

n.kids <- prepped.data[,"n.kids"]
snoob <- prepped.data[,"phi"]
age <- integer(0)
male <- integer(0)

census.data <- prepped.data[,"census"]
census.list <- unique(prepped.data[,"census"])

for(i in 1:length(census.list)){
	census.names <- prepped.data[prepped.data[,"census"]==census.list[i], "id"]
	raw.data.subset <- raw.data[raw.data[,"census"]==census.list[i],]
	male <- c(male, raw.data.subset[match(census.names, raw.data.subset[,"id"]),"male"])
	age <- c(age, raw.data.subset[match(census.names, raw.data.subset[,"id"]),"age"])
}


# drop <- which(census.data==53)
# n.kids <- n.kids[-drop]
# snoob <- snoob[-drop]
# age <- age[-drop]
# male <- male[-drop]
# census.data <- census.data[-drop]

census.n <- tapply(snoob, census.data, length)
current.n <- census.n[as.character(census.data)]
census.phibar <- tapply(snoob, census.data, mean)
current.phibar <- census.phibar[as.character(census.data)]
age.yr <- age/365


# restrict data to childbearing years...will you be between 15 and 60 in the next 5 years?

drop <- which(age < 11*365 | age > 60*365)

n.kids <- n.kids[-drop]
snoob <- snoob[-drop]
age <- age[-drop]
male <- male[-drop]
census.data <- census.data[-drop]

current.n <- current.n[-drop]
current.phibar <- current.phibar[-drop]

age.yr <- age/365


m0p <- mle2(y ~ dpois(lambda=b0), data=list(y=n.kids), start=list(b0=0.2), trace=FALSE)
m0z <- mle2(y ~ my.dzpois(lambda=b0, a=alpha), data=list(y=n.kids), start=list(b0=0.2, alpha=0.1), trace=FALSE)
m0g <- mle2(y ~ dgeom(prob=b0), data=list(y=n.kids), start=list(b0=0.2), trace=FALSE)

m1p <- mle2(y ~ dpois(lambda=exp(b0 + b1*cen(age.yr))), data=list(y=n.kids), start=list(b0=0.2, b1=0), trace=FALSE)
m1z <- mle2(y ~ my.dzpois(lambda=exp(b0 + b1*cen(age.yr)), a=alpha), data=list(y=n.kids), start=list(b0=0.2, b1=0, alpha=0), trace=FALSE)
m1g <- mle2(y ~ dgeom(prob=1/(1+exp(b0 + b1*cen(age.yr)))), data=list(y=n.kids), start=list(b0=0, b1=0), trace=FALSE)

m2p <- mle2(y ~ dpois(lambda=exp(b0 + b1*cen(age.yr) + b2*(cen(age.yr))^2)), data=list(y=n.kids), start=list(b0=0.2, b1=0, b2=0), trace=FALSE)
m2z <- mle2(y ~ my.dzpois(lambda=exp(b0 + b1*cen(age.yr) + b2*(cen(age.yr))^2), a=alpha), data=list(y=n.kids), start=list(b0=0.2, b1=0, b2=0, alpha=0.1), trace=FALSE)
m2g <- mle2(y ~ dgeom(prob=1/(1+exp(b0 + b1*cen(age.yr) + b2*(cen(age.yr))^2))), data=list(y=n.kids), start=list(b0=0.2, b1=0, b2=0), trace=FALSE)

m3p <- mle2(y ~ dpois(lambda=exp(b0 + b1*cen(age.yr) + b2*(cen(age.yr))^2 + b3*snoob)), data=list(y=n.kids), start=list(b0=0.2, b1=0.01, b2=0, b3=0), trace=FALSE)
m3z <- mle2(y ~ my.dzpois(lambda=exp(b0 + b1*cen(age.yr) + b2*(cen(age.yr))^2 + b3*snoob), a=alpha), data=list(y=n.kids), start=list(b0=0.2, b1=0.01, b2=0, b3=0, alpha=0.1), trace=FALSE)
m3g <- mle2(y ~ dgeom(prob=1/(1+exp(b0 + b1*cen(age.yr) + b2*(cen(age.yr))^2 + b3*snoob))), data=list(y=n.kids), start=list(b0=0.0, b1=0.00, b2=0, b3=0), trace=FALSE)

m4p <- mle2(y ~ dpois(lambda=exp(b0 + b1*cen(age.yr) + b2*(cen(age.yr))^2 + b3*snoob + b4*male)), data=list(y=n.kids), start=list(b0=0.2, b1=0.01, b2=0, b3=0, b4=0), trace=FALSE)
m4z <- mle2(y ~ my.dzpois(lambda=exp(b0 + b1*cen(age.yr) + b2*(cen(age.yr))^2 + b3*snoob + b4*male), a=alpha), data=list(y=n.kids), start=list(b0=0.2, b1=0.01, b2=0, b3=0, b4=0, alpha=0.1), trace=FALSE)
m4g <- mle2(y ~ dgeom(prob=1/(1+exp(b0 + b1*cen(age.yr) + b2*(cen(age.yr))^2 + b3*snoob + b4*male))), data=list(y=n.kids), start=list(b0=0.2, b1=0.01, b2=0, b3=0, b4=0), trace=FALSE)

AICtab(m4p, m4z, m4g, m3p, m3z, m3g, m2p, m2z, m2g, m1p, m1z, m1g, m0p, m0z, m0g, weights=T)

# with density dependence

m5p <- mle2(y ~ dpois(lambda=exp((1-current.n/k)*(b0 + b1*cen(age.yr) + b2*(cen(age.yr))^2 + b3*snoob + b4*male))), data=list(y=n.kids), start=list(b0=0.2, b1=0.01, b2=0, b3=0, b4=0, k=10000), trace=FALSE)

m6p <- mle2(y ~ dpois(lambda=exp(exp(-k*current.n)*(b0 + b1*cen(age.yr) + b2*(cen(age.yr))^2 + b3*snoob + b4*male))), data=list(y=n.kids), start=list(b0=0.2, b1=0.01, b2=0, b3=0, b4=0, k=0.001), trace=FALSE)

# with frequency-dependence

m7p <- mle2(y ~ dpois(lambda=exp(b0 + b1*cen(age.yr) + b2*current.phibar)), data=list(y=n.kids), start=list(b0=0.2, b1=0, b2=0), trace=FALSE)
m7z <- mle2(y ~ my.dzpois(lambda=exp(b0 + b1*cen(age.yr) + b2*current.phibar), a=alpha), data=list(y=n.kids), start=list(b0=0.2, b1=0, b2=0, alpha=0.1), trace=FALSE)
m7g <- mle2(y ~ dgeom(prob=1/(1+exp(b0 + b1*cen(age.yr) + b2*current.phibar))), data=list(y=n.kids), start=list(b0=0.2, b1=0, b2=0), trace=FALSE)

m8p <- mle2(y ~ dpois(lambda=exp(b0 + b1*cen(age.yr) + b2*(cen(age.yr))^2 + b3*current.phibar)), data=list(y=n.kids), start=list(b0=0.2, b1=0, b2=0, b3=0), trace=FALSE)
m8z <- mle2(y ~ my.dzpois(lambda=exp(b0 + b1*cen(age.yr) + b2*(cen(age.yr))^2 + b3*current.phibar), a=alpha), data=list(y=n.kids), start=list(b0=0.2, b1=0, b2=0, b3=0, alpha=0.1), trace=FALSE)
m8g <- mle2(y ~ dgeom(prob=1/(1+exp(b0 + b1*cen(age.yr) + b2*(cen(age.yr))^2 + b3*current.phibar))), data=list(y=n.kids), start=list(b0=0.2, b1=0, b2=0, b3=0), trace=FALSE)

m9p <- mle2(y ~ dpois(lambda=exp(b0 + b1*cen(age.yr) + b2*(cen(age.yr))^2 + b3*snoob + b4*current.phibar)), data=list(y=n.kids), start=list(b0=0.2, b1=0.01, b2=0, b3=0, b4=0), trace=FALSE)
m9z <- mle2(y ~ dpois(lambda=exp(b0 + b1*cen(age.yr) + b2*(cen(age.yr))^2 + b3*snoob + b4*current.phibar)), data=list(y=n.kids), start=list(b0=0.2, b1=0.01, b2=0, b3=0, b4=0), trace=FALSE)
m9g <- mle2(y ~ dpois(lambda=exp(b0 + b1*cen(age.yr) + b2*(cen(age.yr))^2 + b3*snoob + b4*current.phibar)), data=list(y=n.kids), start=list(b0=0.2, b1=0.01, b2=0, b3=0, b4=0), trace=FALSE)

m10p <- mle2(y ~ dpois(lambda=exp(b0 + b1*cen(age.yr) + b2*(cen(age.yr))^2 + b3*snoob + b4*male + b5*current.phibar)), data=list(y=n.kids), start=list(b0=0.2, b1=0.01, b2=0, b3=0, b4=0, b5=0), trace=FALSE)
m10z <- mle2(y ~ my.dzpois(lambda=exp(b0 + b1*cen(age.yr) + b2*(cen(age.yr))^2 + b3*snoob + b4*male + b5*current.phibar), a=alpha), data=list(y=n.kids), start=list(b0=-1.15, b1=-0.02, b2=0, b3=0.68, b4=0.03, b5=0.22, alpha=0.1), trace=FALSE)
m10g <- mle2(y ~ dgeom(prob=1/(1+exp(b0 + b1*cen(age.yr) + b2*(cen(age.yr))^2 + b3*snoob + b4*male + b5*current.phibar))), data=list(y=n.kids), start=list(b0=0.2, b1=0.01, b2=0, b3=0, b4=0, b5=0), trace=FALSE)

AICtab(m10p, m10z, m10g, m9p, m9z, m9g, m8p, m8z, m8g, m7p, m7z, m7g, m4p, m4z, m4g, m3p, m3z, m3g, m2p, m2z, m2g, m1p, m1z, m1g, m0p, m0z, m0g, weights=T)

# interactions for the best models, m4 and m10

m4p.i <- mle2(y ~ dpois(lambda=exp(b0 + b1*cen(age.yr) + b2*(cen(age.yr))^2 + b3*snoob + b4*male + b5*male*snoob)), data=list(y=n.kids), start=list(b0=0.2, b1=0.01, b2=0, b3=0, b4=0, b5=0), trace=FALSE)
m4z.i <- mle2(y ~ my.dzpois(lambda=exp(b0 + b1*cen(age.yr) + b2*(cen(age.yr))^2 + b3*snoob + b4*male + b5*male*snoob), a=alpha), data=list(y=n.kids), start=list(b0=0.2, b1=0.01, b2=0, b3=0, b4=0, b5=0, alpha=0.1), trace=FALSE)
m4g.i <- mle2(y ~ dgeom(prob=1/(1+exp(b0 + b1*cen(age.yr) + b2*(cen(age.yr))^2 + b3*snoob + b4*male + b5*male*snoob))), data=list(y=n.kids), start=list(b0=0.2, b1=0.01, b2=0, b3=0, b4=0, b5=0), trace=FALSE)

m10p.i <- mle2(y ~ dpois(lambda=exp(b0 + b1*cen(age.yr) + b2*(cen(age.yr))^2 + b3*snoob + b4*male + b5*current.phibar + b6*male*snoob)), data=list(y=n.kids), start=list(b0=0.2, b1=0.01, b2=0, b3=0, b4=0, b5=0, b6=0), trace=FALSE)
m10z.i <- mle2(y ~ my.dzpois(lambda=exp(b0 + b1*cen(age.yr) + b2*(cen(age.yr))^2 + b3*snoob + b4*male + b5*current.phibar + b6*male*snoob), a=alpha), data=list(y=n.kids), start=list(b0=-1.15, b1=-0.02, b2=0, b3=0.68, b4=0.03, b5=0.22, b6=0, alpha=0.1), trace=FALSE)
m10g.i <- mle2(y ~ dgeom(prob=1/(1+exp(b0 + b1*cen(age.yr) + b2*(cen(age.yr))^2 + b3*snoob + b4*male + b5*current.phibar + b6*male*snoob))), data=list(y=n.kids), start=list(b0=0.2, b1=0.01, b2=0, b3=0, b4=0, b5=0, b6=0), trace=FALSE)

AICtab(m10p.i, m10z.i, m10g.i, m4p.i, m4z.i, m4g.i, m10p, m10z, m10g, m9p, m9z, m9g, m8p, m8z, m8g, m7p, m7z, m7g, m4p, m4z, m4g, m3p, m3z, m3g, m2p, m2z, m2g, m1p, m1z, m1g, m0p, m0z, m0g, weights=T)

# compare(m10p.i, m10z.i, m10g.i, m4p.i, m4z.i, m4g.i, m10p, m10z, m10g, m9p, m9z, m9g, m8p, m8z, m8g, m7p, m7z, m7g, m4p, m4z, m4g, m3p, m3z, m3g, m2p, m2z, m2g, m1p, m1z, m1g, m0p, m0z, m0g, nobs=length(n.kids))








# # now let's graph the winners

b0 <- -0.878
b1 <- -0.006
b2 <- -0.002
b3 <- 0.681
b4 <- 0.120
b5 <- 0.244
alpha <- 0.357

# my.dzpois(lambda=exp(b0 + b1*cen(age.yr) + b2*(cen(age.yr))^2 + b3*snoob + b4*male + b5*current.phibar + b6*male*snoob), a=alpha)



png("./temp/fertility.png", height=15, width=20, 
        units='cm', res=300)

par(mfrow=c(2,2))

k <- 1
plot.phibar <- 0.5
will.have.k <- as.numeric(n.kids==k)
age.yr[age.yr >= 60] <- 59
s.pr.will.have.k <- tapply(will.have.k[snoob==1 & abs(current.phibar - plot.phibar) < 0.1], floor(age.yr[snoob==1 & abs(current.phibar - 0.5) < 0.1]), mean)
ns.pr.will.have.k <- tapply(will.have.k[snoob==0 & abs(current.phibar - plot.phibar) < 0.1], floor(age.yr[snoob==0 & abs(current.phibar - 0.5) < 0.1]), mean)
age.list <- 11:59

plot(ns.pr.will.have.k ~ age.list, ylim=c(min(c(ns.pr.will.have.k, s.pr.will.have.k)), max(ns.pr.will.have.k, s.pr.will.have.k)), ylab="", xlab="Age", las=1, main=paste("Prob. of having", k, "kid", sep=" "), col="gray")
points(s.pr.will.have.k ~ age.list, col="black", pch=1)

curve(my.dzpois(k, lambda=exp(b0 + b1*(x-32.67694) + b2*(x-32.67694)^2 + b3*1 + b4*1 + b5*plot.phibar), a=alpha), from=11, to=60, col="dodgerblue", ylim=c(0, 0.2), add=T)
curve(my.dzpois(k, lambda=exp(b0 + b1*(x-32.67694) + b2*(x-32.67694)^2 + b3*0 + b4*1 + b5*plot.phibar), a=alpha), from=11, to=60, col="dodgerblue", lty=1, add=T)
curve(my.dzpois(k, lambda=exp(b0 + b1*(x-32.67694) + b2*(x-32.67694)^2 + b3*1 + b4*0 + b5*plot.phibar), a=alpha), from=11, to=60, col="maroon1", add=T)
curve(my.dzpois(k, lambda=exp(b0 + b1*(x-32.67694) + b2*(x-32.67694)^2 + b3*0 + b4*0 + b5*plot.phibar), a=alpha), from=11, to=60, col="maroon1", lty=1, add=T)


k <- 2
plot.phibar <- 0.5
will.have.k <- as.numeric(n.kids==k)
s.pr.will.have.k <- tapply(will.have.k[snoob==1 & abs(current.phibar - plot.phibar) < 0.1], floor(age.yr[snoob==1 & abs(current.phibar - 0.5) < 0.1]), mean)
ns.pr.will.have.k <- tapply(will.have.k[snoob==0 & abs(current.phibar - plot.phibar) < 0.1], floor(age.yr[snoob==0 & abs(current.phibar - 0.5) < 0.1]), mean)
age.list <- 11:59

plot(ns.pr.will.have.k ~ age.list, ylim=c(min(c(ns.pr.will.have.k, s.pr.will.have.k)), max(ns.pr.will.have.k, s.pr.will.have.k)), ylab="", xlab="Age", las=1, main=paste("Prob. of having", k, "kids", sep=" "), col="gray")
points(s.pr.will.have.k ~ age.list, col="black", pch=1)

curve(my.dzpois(k, lambda=exp(b0 + b1*(x-32.67694) + b2*(x-32.67694)^2 + b3*1 + b4*1 + b5*plot.phibar), a=alpha), from=11, to=60, col="dodgerblue", ylim=c(0, 0.2), add=T)
curve(my.dzpois(k, lambda=exp(b0 + b1*(x-32.67694) + b2*(x-32.67694)^2 + b3*0 + b4*1 + b5*plot.phibar), a=alpha), from=11, to=60, col="dodgerblue", lty=1, add=T)
curve(my.dzpois(k, lambda=exp(b0 + b1*(x-32.67694) + b2*(x-32.67694)^2 + b3*1 + b4*0 + b5*plot.phibar), a=alpha), from=11, to=60, col="maroon1", add=T)
curve(my.dzpois(k, lambda=exp(b0 + b1*(x-32.67694) + b2*(x-32.67694)^2 + b3*0 + b4*0 + b5*plot.phibar), a=alpha), from=11, to=60, col="maroon1", lty=1, add=T)


k <- 3
plot.phibar <- 0.5
will.have.k <- as.numeric(n.kids==k)
s.pr.will.have.k <- tapply(will.have.k[snoob==1 & abs(current.phibar - plot.phibar) < 0.1], floor(age.yr[snoob==1 & abs(current.phibar - 0.5) < 0.1]), mean)
ns.pr.will.have.k <- tapply(will.have.k[snoob==0 & abs(current.phibar - plot.phibar) < 0.1], floor(age.yr[snoob==0 & abs(current.phibar - 0.5) < 0.1]), mean)
age.list <- 11:59

plot(ns.pr.will.have.k ~ age.list, ylim=c(min(c(ns.pr.will.have.k, s.pr.will.have.k)), max(ns.pr.will.have.k, s.pr.will.have.k)), ylab="", xlab="Age", las=1, main=paste("Prob. of having", k, "kids", sep=" "), col="gray")
points(s.pr.will.have.k ~ age.list, col="black", pch=1)

curve(my.dzpois(k, lambda=exp(b0 + b1*(x-32.67694) + b2*(x-32.67694)^2 + b3*1 + b4*1 + b5*plot.phibar), a=alpha), from=11, to=60, col="dodgerblue", ylim=c(0, 0.2), add=T)
curve(my.dzpois(k, lambda=exp(b0 + b1*(x-32.67694) + b2*(x-32.67694)^2 + b3*0 + b4*1 + b5*plot.phibar), a=alpha), from=11, to=60, col="dodgerblue", lty=1, add=T)
curve(my.dzpois(k, lambda=exp(b0 + b1*(x-32.67694) + b2*(x-32.67694)^2 + b3*1 + b4*0 + b5*plot.phibar), a=alpha), from=11, to=60, col="maroon1", add=T)
curve(my.dzpois(k, lambda=exp(b0 + b1*(x-32.67694) + b2*(x-32.67694)^2 + b3*0 + b4*0 + b5*plot.phibar), a=alpha), from=11, to=60, col="maroon1", lty=1, add=T)


k <- 4
plot.phibar <- 0.5
will.have.k <- as.numeric(n.kids==k)
s.pr.will.have.k <- tapply(will.have.k[snoob==1 & abs(current.phibar - plot.phibar) < 0.1], floor(age.yr[snoob==1 & abs(current.phibar - 0.5) < 0.1]), mean)
ns.pr.will.have.k <- tapply(will.have.k[snoob==0 & abs(current.phibar - plot.phibar) < 0.1], floor(age.yr[snoob==0 & abs(current.phibar - 0.5) < 0.1]), mean)
age.list <- 11:59

plot(ns.pr.will.have.k ~ age.list, ylim=c(min(c(ns.pr.will.have.k, s.pr.will.have.k)), max(ns.pr.will.have.k, s.pr.will.have.k)), ylab="", xlab="Age", las=1, main=paste("Prob. of having", k, "kids", sep=" "), col="gray")
points(s.pr.will.have.k ~ age.list, col="black", pch=1)

curve(my.dzpois(k, lambda=exp(b0 + b1*(x-32.67694) + b2*(x-32.67694)^2 + b3*1 + b4*1 + b5*plot.phibar), a=alpha), from=11, to=60, col="dodgerblue", ylim=c(0, 0.2), add=T)
curve(my.dzpois(k, lambda=exp(b0 + b1*(x-32.67694) + b2*(x-32.67694)^2 + b3*0 + b4*1 + b5*plot.phibar), a=alpha), from=11, to=60, col="dodgerblue", lty=1, add=T)
curve(my.dzpois(k, lambda=exp(b0 + b1*(x-32.67694) + b2*(x-32.67694)^2 + b3*1 + b4*0 + b5*plot.phibar), a=alpha), from=11, to=60, col="maroon1", add=T)
curve(my.dzpois(k, lambda=exp(b0 + b1*(x-32.67694) + b2*(x-32.67694)^2 + b3*0 + b4*0 + b5*plot.phibar), a=alpha), from=11, to=60, col="maroon1", lty=1, add=T)

dev.off()














# Individual Change

# outcome: will you be a snoob next census?

drop <- which(prepped.data[,"state"]==3)
living.data <- prepped.data[-drop,]

census.list <- unique(living.data[,"census"])
census.data <- living.data[,"census"]

is.snoob <- living.data[,"phi"]

will.be.snoob <- rep(0, nrow(living.data))
age <- integer(0)
male <- integer(0)

for(i in 1:(length(census.list)-1)){
	census.names <- living.data[living.data[,"census"]==census.list[i], "id"]
	prepped.subset <- prepped.data[prepped.data[,"census"]==census.list[i+1],]
	will.be.snoob[which(census.data==census.list[i])] <- prepped.subset[match(census.names, prepped.subset[,"id"]),"phi"]
	
	raw.data.subset <- raw.data[raw.data[,"census"]==census.list[i],]
	male <- c(male, raw.data.subset[match(census.names, raw.data.subset[,"id"]),"male"])
	age <- c(age, raw.data.subset[match(census.names, raw.data.subset[,"id"]),"age"])
}

phibars <- tapply(prepped.data[,"phi"], prepped.data[,"census"], mean)
census.phibar <- phibars[as.character(living.data[,"census"])]

drop <- which(census.data==53)
is.snoob <- is.snoob[-drop]
will.be.snoob <- will.be.snoob[-drop]
census.phibar <- census.phibar[-drop]

pop.size <- tapply(prepped.data[,"phi"], prepped.data[,"census"], length)
current.popsize <- pop.size[as.character(living.data[,"census"])]
current.popsize <- current.popsize[-drop]
current.census <- census.data[-drop]







age.yr <- age/365

m0 <- mle2(y ~ dbinom(1, prob=p), data=list(y=is.snoob), start=list(p=0.2), trace=FALSE)

m1 <- mle2(y ~ dbinom(1, prob=L*(x^(B)/(x^(B) + (1-x)^(B))) + (1-L)*z), data=list(y = will.be.snoob, x =census.phibar, z=is.snoob), start=list(L=0.2, B=0.2), trace=FALSE)

m2 <- mle2(y ~ dbinom(1, prob=L*(x + 4*(B-1)*x*(1-x)*(x-0.5))+ (1-L)*z), data=list(y = will.be.snoob, x =census.phibar, z=is.snoob), start=list(L=0.2, B=0), trace=FALSE)

m3 <- mle2(y ~ dbinom(1, prob=L*(x + 4*(B-1)*x*(1-x)*(x-k))+ (1-L)*z), data=list(y = will.be.snoob, x =census.phibar, z=is.snoob), start=list(L=0.2, B=0, k=0.5), trace=FALSE)

m4 <-  mle2(y ~ dbinom(1, prob=L*(exp(b0 + b1*x)/(1+exp(b0 + b1*x))) + (1-L)*z ), data=list(y = will.be.snoob, x = census.phibar, z = is.snoob), start=list(L=0.2, b0=0.3, b1=0), trace=FALSE)

m5 <- mle2(y ~ dbinom(1, prob=L*(1/(1+exp(b0 + b1*n))) + (1-L)*z), data=list(y = will.be.snoob, n = current.popsize, z = is.snoob), start=list(L=0.2, b0=0.3, b1=0), trace=FALSE)

m6 <- mle2(y ~ dbinom(1, prob=L*(1/(1+exp(b0 + b1*su(n) + b2*su(n)^2))) + (1-L)*z), data=list(y = will.be.snoob, n = current.popsize, z = is.snoob), start=list(L=0.2, b0=0.3, b1=0, b2=0), trace=FALSE)

m7 <- mle2(y ~ dbinom(1, prob=L*(1/(1+exp(b0 + b1*t))) + (1-L)*z ), data=list(y = will.be.snoob, t = current.census, z = is.snoob), start=list(L=0.2, b0=0.3, b1=0), trace=FALSE)

m8 <- mle2(y ~ dbinom(1, prob=L*(1/(1+exp(b0 + b1*t + b2*t^2))) + (1-L)*z ), data=list(y = will.be.snoob, t = current.census, z = is.snoob), start=list(L=0.2, b0=0.3, b1=0, b2=0), trace=FALSE)

m9 <- mle2(y ~ dbinom(1, prob= L*(b0 + b1*sin(b2*t)) + (1-L)*z), data=list(y = will.be.snoob, t = current.census, z = is.snoob), start=list(L=0.2, b0=0.3, b1=0.1, b2=1), trace=FALSE)

m10 <- mle2(y ~ dbinom(1, prob= L*(b0 + b1*sin(b2*n)) + (1-L)*z), data=list(y = will.be.snoob, n = current.popsize, z = is.snoob), start=list(L=0.2, b0=0.3, b1=0.1, b2=1), trace=FALSE)

m11 <- mle2(y ~ dbinom(1, prob=L*(1/(1+exp(b0 + b1*n + b2*t))) + (1-L)*z), data=list(y = will.be.snoob, n = current.popsize, t = current.census, z = is.snoob), start=list(L=0.2, b0=0.3, b1=0, b2=0), trace=FALSE)

m12 <-  mle2(y ~ dbinom(1, prob=L*(1/(1+exp(b0 + b1*x + b2*a + b3*m))) + (1-L)*z ), data=list(y = will.be.snoob, x = census.phibar, z = is.snoob, a = cen(age.yr), m = male), start=list(L=0.2, b0=0.3, b1=0, b2=0, b3=0), trace=FALSE)

m13 <-  mle2(y ~ dbinom(1, prob=L*(1/(1+exp(b0 + b1*x + b2*a + b3*m + b4*t))) + (1-L)*z ), data=list(y = will.be.snoob, x = census.phibar, a = cen(age.yr), m = male, z = is.snoob, t = current.census), start=list(L=0.2, b0=0.3, b1=0, b2=0, b3=0, b4=0), trace=FALSE)
 
# compare(m0,m1,m2,m3,m4,m5,m6,m7,m8, m9, m10, m11, m12,m13, nobs=length(will.be.snoob))







census.tag <- census.data[census.data!=0]

bin1 <- which(is.snoob==0)
bin2 <- which(is.snoob==1)

non.to.snoob <- tapply(will.be.snoob[bin1], census.tag[bin1], mean)
snoob.to.snoob <- tapply(will.be.snoob[bin2], census.tag[bin2], mean)








png("./temp/learning.png", height=15, width=20, 
        units='cm', res=300)

par(mfrow=c(2,1))

plot(snoob.to.snoob ~ phibars[-1], ylim=c(0.85, 1), ylab="Prob. Remaining Snoob", xlab="Mean Pop. Frq", col="blue", las=1, frame.plot=F, xlim=c(.15, 1))

# L * ( x^(B) / ( x^(B) + (1-x)^(B)) ) + (1-L) * z ) 
# data=list(y = will.be.snoob, x = census.phibar, z=is.snoob)

# L <- 0.104
# D <- 0.619
# curve( L * ( x + D*x*(1-x)*(2*x-1) )+ (1-L), from=0, to=1, col="dodgerblue", add=T)

L <- 0.09996
D <- 1.37388
curve( L * ( x^D / (x^D + (1-x)^D) ) + ( 1 - L ), from=0, to=1, ylim=c(0,1), col="black", add=T)	


post <- sample.naive.posterior(m1)
frqs <- seq(0.01, 0.99, length=100)

s.p.lb <- rep(NA, length(frqs))
s.p.ub <- rep(NA, length(frqs))

# for(i in 1:length(frqs)){
# 	prob <- post$L*(frqs[i])^(post$B)/(frqs[i]^post$B + (1-frqs[i])^post$B)
# 	s.p.lb[i] <- HPDI(prob)[1]
# 	s.p.ub[i] <- HPDI(prob)[2]
# }

# points(frqs, s.p.lb, lty=1, type="l")
# points(frqs, s.p.ub, lty=1, type="l")



plot(non.to.snoob ~ phibars[-1], ylim=c(0,.15), ylab="Prob. Becoming Snoob", xlab="Mean Pop. Frq", col="blue", las=1, frame.plot=F, xlim=c(.15, 1))

# L <- 0.104
# D <- 0.619
# curve(L*(x + D*x*(1-x)*(2*x-1)), from=0, to=1, ylim=c(0, 1), col="dodgerblue", add=T)

L <- 0.09996
D <- 1.37388
curve(L*(x^D/(x^D + (1-x)^D)), from=0, to=1, ylim=c(0,1), col="black", add=T)


# post <- sample.naive.posterior(m1)
# frqs <- seq(0, 1, length=100)

# ns.p.ci <- sapply(frqs, function(z) HPDI( post$L*(z)^(post$B)/(z^post$B + (1-z)^post$B)))
# points(frqs, ns.p.ci[1,], lty=1, type="l")
# points(frqs, ns.p.ci[2,], lty=1, type="l")



dev.off()

dir_init('./output')

files <- list.files('./temp', full.names=TRUE)
file.copy(files, './output')

print(paste(Sys.time(), "- analysis complete"))