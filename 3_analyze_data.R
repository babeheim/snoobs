
rm(list = ls())

source("project_support.R")

d <- read.csv("analysis_data.csv")



cat("prep predictors\n")

# for some reason we call the initial census 'census 0' which is weird
census_list <- unique(d$census)
n_censuses <- length(census_list)
stopifnot(census_list == (seq_len(n_censuses) - 1))

cens <- data.frame(
  census = census_list
)

# outcome: will you be a snoob next census?

d <- rename(d, snoob = phi, n_kids = n.kids)

# state 3 in the prepped data means 'died between the last census and the current one'
d$died <- as.numeric(d$state == 3)
alive <- which(d$died != 1)

cens$phi_bar <- as.numeric(tapply(d$snoob[alive], d$census[alive], mean))
d$phi_bar <- cens$phi_bar[match(d$census, cens$census)]

cens$pop_size <- as.numeric(table(d$census[alive]))
d$pop_size <- cens$pop_size[match(d$census, cens$census)]

d$will_die <- NA
d$will_be_snoob <- NA

for (i in 1:(n_censuses - 1)) {
  this_census <- which(d$census == cens$census[i])
  died_next_census <- which(d$died == 1 & d$census == cens$census[i + 1])
  d$will_die[this_census] <- as.numeric(d$id[this_census] %in% d$id[died_next_census])
  snoob_next_census <- which(d$snoob == 1 & d$census == cens$census[i + 1])
  d$will_be_snoob[this_census] <- as.numeric(d$id[this_census] %in% d$id[snoob_next_census])
}

nonsnoobs <- which(d$snoob == 0)
snoobs <- which(d$snoob == 1)

cens$non_to_snoob <- tapply(d$will_be_snoob[nonsnoobs], d$census[nonsnoobs], mean)
cens$snoob_to_snoob <- tapply(d$will_be_snoob[snoobs], d$census[snoobs], mean)

# exclude final census, because we don't know what will happen after it
keep <- which(d$census != census_list[n_censuses])
d <- d[keep,]

stopifnot(!any(is.na(d$will_die)))
stopifnot(!any(is.na(d$age)))
stopifnot(!any(is.na(d$male)))

d$age_c <- floor(d$age)
d$age_c[d$age_c > 100] <- 100



cat("fit mortality models\n")

cat("fit mortality model m0\n")
m0 <- mle2(minuslogl = y ~ dbinom(1, prob=1/(1+exp(a))),
  data = list(y = d$will_die), start = list(a = 0.1), trace = FALSE)

cat("fit mortality model m1\n")
m1 <- mle2(minuslogl = y ~ dbinom(1, prob=1/(1+exp(p1 + p2*cen(d$age)))),
  data = list(y = d$will_die), start=list(p1 = 4.8, p2 = 0), trace = FALSE)

cat("fit mortality model m2\n")
m2 <- mle2(minuslogl = y ~ dbinom(1, prob=1/(1+exp(p1 + p2*cen(d$age) + p3*cen(d$age)^2))),
  data = list(y = d$will_die), start = list(p1 = 4.8, p2 = 0, p3 = 0), trace = FALSE)

cat("fit mortality model m3\n")
m3 <- mle2(y ~ dbinom(1, prob=1/(1+exp(p1 + p2*cen(d$age) + p3*cen(d$age)^2 + p4*d$male))),
  data=list(y = d$will_die), start=list(p1=4.8, p2=0, p3=0, p4=0), trace=FALSE)

cat("fit mortality model m4\n")
m4 <- mle2(y ~ dbinom(1, prob=1/(1+exp(p1 + p2*cen(d$age) + p3*cen(d$age)^2 + p4*d$male + p5*d$snoob))),
  data=list(y = d$will_die), start=list(p1=4.8, p2=0, p3=0, p4=0, p5=0), trace=FALSE) # warnings

cat("fit mortality model m5\n")
m5 <- mle2(
  minuslogl = y ~ dbinom(1, prob=1/(1+exp(p1 + p2*cen(d$age) + p3*cen(d$age)^2 + p4*d$male + p5*d$snoob + p6*d$phi_bar))),
  data = list(y = d$will_die), start = list(p1=4.8, p2=0, p3=0, p4=0, p5=0, p6=0), trace = FALSE
) # warnings

cat("fit mortality model m6\n")
m6 <- mle2(y ~ dbinom(1, prob=1/(1+exp(p1 + p2*cen(d$age) + p3*cen(d$age)^2 + p4*d$male + p5*d$snoob + p6*su(d$age)^3))),
  data=list(y = d$will_die), start=list(p1=4.8, p2=0, p3=0, p4=0, p5=0, p6=0), trace=FALSE) # warnings

cat("fit mortality model m0z\n")
m0z <- mle2(y ~ my.dzbinom(1, prob=1/(1+exp(p1)), a=alpha),
  data=list(y = d$will_die), start=list(p1=0.1, alpha=0.1), trace=FALSE) # warnings

cat("fit mortality model m1z\n")
m1z <- mle2(y ~ my.dzbinom(1, prob=1/(1+exp(p1 + p2*cen(d$age))), a=alpha),
  data=list(y = d$will_die), start=list(p1=3.14, p2=0, alpha=0.1), trace=FALSE) # warnings

# cat("fit mortality model m2z\n")
# m2z <- mle2(y ~ my.dzbinom(1, prob=1/(1+exp(p1 + p2*cen(d$age) + p3*cen(d$age)^2)), a=alpha) , data=list(y = d$will_die), start=list(p1=3.66, p2=0, p3=0, alpha=0.1), trace=FALSE)

cat("fit mortality model m3z\n")
m3z <- mle2(y ~ my.dzbinom(1, prob=1/(1+exp(p1 + p2*cen(d$age) + p3*cen(d$age)^2 + p4*d$male)), a=alpha),
  data=list(y = d$will_die), start=list(p1=5.39, p2=0, p3=0, p4=0, alpha=0.1), trace=FALSE) # warnings

cat("fit mortality model m4z\n")
m4z <- mle2(y ~ my.dzbinom(1, prob=1/(1+exp(p1 + p2*cen(d$age) + p3*cen(d$age)^2 + p4*d$male + p5*d$snoob)), a=1/(1+exp(alpha))),
  data=list(y = d$will_die), start=list(p1=4.8, p2=0, p3=0, p4=0, p5=0, alpha=0.1), trace=FALSE) # warnings

cat("fit mortality model m5z\n")
m5z <- mle2(y ~ my.dzbinom(1, prob=1/(1+exp(p1 + p2*cen(d$age) + p3*cen(d$age)^2 + p4*d$male + p5*d$snoob + p6*d$phi_bar)), a=1/(1+exp(alpha))),
  data=list(y = d$will_die), start=list(p1=4.8, p2=0, p3=0, p4=0, p5=0, p6=0, alpha=0.1), trace=FALSE) # warnings

AICtab(m0, m1, m2, m3, m4, m5, m0z, m1z, m3z, m4z, m5z, m6, weights = TRUE)

cat("create mortality.png\n")

png("./figures/mortality.png", height=15, width=20, 
        units='cm', res=300)

par(mfrow=c(2,2))

s.pr.die.nums <- tapply(d$will_die[d$snoob==1], d$age_c[d$snoob==1], mean)
s.pr.die <- rep(0, length(0:100))
names(s.pr.die) <- 0:100
s.pr.die[names(s.pr.die.nums)] <- s.pr.die.nums

ns.pr.die.nums <- tapply(d$will_die[d$snoob==0], d$age_c[d$snoob==0], mean)
ns.pr.die <- rep(0, length(0:100))
names(ns.pr.die) <- 0:100
ns.pr.die[names(ns.pr.die.nums)] <- ns.pr.die.nums

# what model is this from??
p1 <- 4.673
p2 <- -0.036
p3 <- -0.001
p4 <- -0.088
p5 <- 0.559

age_cen <- 33.81

plot(1000*ns.pr.die ~ names(ns.pr.die), las=1)
points(1000*s.pr.die ~ names(ns.pr.die), pch=20)
curve(1000*dbinom(1, 1, prob=1/(1+exp(p1 + p2*(x-age_cen) + p3*(x-age_cen)^2 + p4*0 + p5*0))), type="l", lty=2, from=0, to=100, add=T, col="maroon1")
curve(1000*dbinom(1, 1, prob=1/(1+exp(p1 + p2*(x-age_cen) + p3*(x-age_cen)^2 + p4*0 + p5*1))), type="l", from=0, to=100, add=T, col="maroon1")
curve(1000*dbinom(1, 1, prob=1/(1+exp(p1 + p2*(x-age_cen) + p3*(x-age_cen)^2 + p4*1 + p5*0))), type="l", lty=2, from=0, to=100, add=T, col="dodgerblue")
curve(1000*dbinom(1, 1, prob=1/(1+exp(p1 + p2*(x-age_cen) + p3*(x-age_cen)^2 + p4*1 + p5*1))), type="l", from=0, to=100, add=T, col="dodgerblue")

s.pr.die.nums <- tapply(d$will_die[d$snoob==1], d$age_c[d$snoob==1], mean)
s.pr.die <- rep(0, length(0:100))
names(s.pr.die) <- 0:100
s.pr.die[names(s.pr.die.nums)] <- s.pr.die.nums

ns.pr.die.nums <- tapply(d$will_die[d$snoob==0], d$age_c[d$snoob==0], mean)
ns.pr.die <- rep(0, length(0:100))
names(ns.pr.die) <- 0:100
ns.pr.die[names(ns.pr.die.nums)] <- ns.pr.die.nums

# from what model?
p1 <- 4.673
p2 <- -0.036
p3 <- -0.001
p4 <- -0.088
p5 <- 0.559

plot(1000*ns.pr.die ~ names(ns.pr.die), log="y", las=1, col="gray", ylim=c(.1, 1000))
points(1000*s.pr.die ~ names(ns.pr.die), pch=1)

x <- c(0:100)
points(1000*dbinom(1, 1, prob=1/(1+exp(p1 + p2*(x-age_cen) + p3*(x-age_cen)^2 + p4*0 + p5*0))), type="l", col="maroon1")
points(1000*dbinom(1, 1, prob=1/(1+exp(p1 + p2*(x-age_cen) + p3*(x-age_cen)^2 + p4*0 + p5*1))), type="l", col="maroon1")
points(1000*dbinom(1, 1, prob=1/(1+exp(p1 + p2*(x-age_cen) + p3*(x-age_cen)^2 + p4*1 + p5*0))), type="l", col="dodgerblue")
points(1000*dbinom(1, 1, prob=1/(1+exp(p1 + p2*(x-age_cen) + p3*(x-age_cen)^2 + p4*1 + p5*1))), type="l", col="dodgerblue")

human.ages <- 0:100
plot(cumprod(1-ns.pr.die))
points(cumprod(1-s.pr.die), pch=20)
fitted.s <- dbinom(1, 1, prob=1/(1+exp(p1 + p2*(human.ages-age_cen) + p3*(human.ages-age_cen)^2 + p4*0 + p5*1)))
fitted.ns <- dbinom(1, 1, prob=1/(1+exp(p1 + p2*(human.ages-age_cen) + p3*(human.ages-age_cen)^2 + p4*0 + p5*0)))
points(cumprod(1-fitted.s), type="l", col="maroon1")
points(cumprod(1-fitted.ns), type="l", col="maroon1")
fitted.s <- dbinom(1, 1, prob=1/(1+exp(p1 + p2*(human.ages-age_cen) + p3*(human.ages-age_cen)^2 + p4*1 + p5*1)))
fitted.ns <- dbinom(1, 1, prob=1/(1+exp(p1 + p2*(human.ages-age_cen) + p3*(human.ages-age_cen)^2 + p4*1 + p5*0)))
points(cumprod(1-fitted.s), type="l", col="dodgerblue")
points(cumprod(1-fitted.ns), type="l", col="dodgerblue")

# Snoobs tend to be younger! 

ns.age <- tapply(d$age[d$snoob==0], d$census[d$snoob==0], mean)
s.age <- tapply(d$age[d$snoob==1], d$census[d$snoob==1], mean)
plot(ns.age/365, type="l", ylim=c(25, 45), ylab="Age in Years", xlab="Census")
points(s.age/365, col="red", type="l")

dev.off()



cat("fit fertility models\n")

# restrict data to childbearing years...will you be between 15 and 60 in the next 5 years?
keep <- which(d$age >= 11 & d$age <= 60)
dk <- d[keep,]

cat("fit fertility model m0p\n")
m0p <- mle2(y ~ dpois(lambda = b0),
  data = list(y = dk$n_kids), start = list(b0 = 0.2), trace = FALSE)

cat("fit fertility model m0z\n")
m0z <- mle2(y ~ my.dzpois(lambda = b0, a = alpha),
  data = list(y = dk$n_kids), start = list(b0 = 0.2, alpha = 0.1), trace = FALSE)

cat("fit fertility model m0g\n")
m0g <- mle2(y ~ dgeom(prob = b0),
  data=list(y = dk$n_kids), start = list(b0 = 0.2), trace = FALSE)

cat("fit fertility model m1p\n")
m1p <- mle2(y ~ dpois(lambda=exp(b0 + b1*cen(dk$age))), data=list(y = dk$n_kids), start=list(b0=0.2, b1=0), trace=FALSE)

cat("fit fertility model m1z\n")
m1z <- mle2(y ~ my.dzpois(lambda = exp(b0 + b1*cen(dk$age)), a = alpha),
  data=list(y = dk$n_kids), start = list(b0 = 0.2, b1 = 0, alpha = 0), trace = FALSE)

cat("fit fertility model m1g\n")
m1g <- mle2(y ~ dgeom(prob = 1/(1 + exp(b0 + b1*cen(dk$age)))),
  data = list(y = dk$n_kids), start = list(b0 = 0, b1 = 0), trace = FALSE)

cat("fit fertility model m2p\n")
m2p <- mle2(y ~ dpois(lambda = exp(b0 + b1*cen(dk$age) + b2*(cen(dk$age))^2)),
  data = list(y = dk$n_kids), start = list(b0 = 0.2, b1 = 0, b2 = 0), trace = FALSE)

cat("fit fertility model m2z\n")
m2z <- mle2(y ~ my.dzpois(lambda = exp(b0 + b1*cen(dk$age) + b2*(cen(dk$age))^2), a = alpha),
  data = list(y = dk$n_kids), start = list(b0 = 0.2, b1 = 0, b2 = 0, alpha = 0.1), trace = FALSE)

cat("fit fertility model m2g\n")
m2g <- mle2(y ~ dgeom(prob=1/(1+exp(b0 + b1*cen(dk$age) + b2*(cen(dk$age))^2))), data=list(y = dk$n_kids), start=list(b0=0.2, b1=0, b2=0), trace=FALSE)

cat("fit fertility model m3p\n")
m3p <- mle2(y ~ dpois(lambda = exp(b0 + b1*cen(dk$age) + b2*(cen(dk$age))^2 + b3*dk$snoob)),
  data = list(y = dk$n_kids), start = list(b0 = 0.2, b1 = 0.01, b2 = 0, b3 = 0), trace = FALSE)

cat("fit fertility model m3z\n")
m3z <- mle2(y ~ my.dzpois(lambda = exp(b0 + b1*cen(dk$age) + b2*(cen(dk$age))^2 + b3*dk$snoob), a=alpha),
  data = list(y = dk$n_kids), start = list(b0 = 0.2, b1 = 0.01, b2 = 0, b3 = 0, alpha = 0.1), trace = FALSE)

cat("fit fertility model m3g\n")
m3g <- mle2(y ~ dgeom(prob=1/(1+exp(b0 + b1*cen(dk$age) + b2*(cen(dk$age))^2 + b3*dk$snoob))),
  data = list(y = dk$n_kids), start = list(b0 = 0.0, b1 = 0.00, b2 = 0, b3 = 0), trace = FALSE)

cat("fit fertility model m4p\n")
m4p <- mle2(y ~ dpois(lambda=exp(b0 + b1*cen(dk$age) + b2*(cen(dk$age))^2 + b3*dk$snoob + b4*dk$male)),
  data = list(y = dk$n_kids), start = list(b0 = 0.2, b1 = 0.01, b2 = 0, b3 = 0, b4 = 0), trace = FALSE)

cat("fit fertility model m4z\n")
m4z <- mle2(y ~ my.dzpois(lambda=exp(b0 + b1*cen(dk$age) + b2*(cen(dk$age))^2 + b3*dk$snoob + b4*dk$male), a=alpha), data=list(y = dk$n_kids), start=list(b0=0.2, b1=0.01, b2=0, b3=0, b4=0, alpha=0.1), trace=FALSE)

cat("fit fertility model m4g\n")
m4g <- mle2(y ~ dgeom(prob=1/(1+exp(b0 + b1*cen(dk$age) + b2*(cen(dk$age))^2 + b3*dk$snoob + b4*dk$male))), data=list(y = dk$n_kids), start=list(b0=0.2, b1=0.01, b2=0, b3=0, b4=0), trace=FALSE)

# AICtab(m4p, m4z, m4g, m3p, m3z, m3g, m2p, m2z, m2g, m1p, m1z, m1g, m0p, m0z, m0g, weights=T)

# with density dependence

cat("fit fertility model m5p\n")
m5p <- mle2(y ~ dpois(lambda=exp((1-dk$pop_size/k)*(b0 + b1*cen(dk$age) + b2*(cen(dk$age))^2 + b3*dk$snoob + b4*dk$male))), data=list(y = dk$n_kids), start=list(b0=0.2, b1=0.01, b2=0, b3=0, b4=0, k=10000), trace=FALSE)

cat("fit fertility model m6p\n")
m6p <- mle2(y ~ dpois(lambda=exp(exp(-k*dk$pop_size)*(b0 + b1*cen(dk$age) + b2*(cen(dk$age))^2 + b3*dk$snoob + b4*dk$male))), data=list(y = dk$n_kids), start=list(b0=0.2, b1=0.01, b2=0, b3=0, b4=0, k=0.001), trace=FALSE)

# with frequency-dependence

cat("fit fertility model m7p\n")
m7p <- mle2(y ~ dpois(lambda=exp(b0 + b1*cen(dk$age) + b2*dk$phi_bar)), data=list(y = dk$n_kids), start=list(b0=0.2, b1=0, b2=0), trace=FALSE)

cat("fit fertility model m7z\n")
m7z <- mle2(y ~ my.dzpois(lambda=exp(b0 + b1*cen(dk$age) + b2*dk$phi_bar), a=alpha), data=list(y = dk$n_kids), start=list(b0=0.2, b1=0, b2=0, alpha=0.1), trace=FALSE)

cat("fit fertility model m7g\n")
m7g <- mle2(y ~ dgeom(prob=1/(1+exp(b0 + b1*cen(dk$age) + b2*dk$phi_bar))), data=list(y = dk$n_kids), start=list(b0=0.2, b1=0, b2=0), trace=FALSE)

cat("fit fertility model m8p\n")
m8p <- mle2(y ~ dpois(lambda=exp(b0 + b1*cen(dk$age) + b2*(cen(dk$age))^2 + b3*dk$phi_bar)), data=list(y = dk$n_kids), start=list(b0=0.2, b1=0, b2=0, b3=0), trace=FALSE)

cat("fit fertility model m8z\n")
m8z <- mle2(y ~ my.dzpois(lambda=exp(b0 + b1*cen(dk$age) + b2*(cen(dk$age))^2 + b3*dk$phi_bar), a=alpha), data=list(y = dk$n_kids), start=list(b0=0.2, b1=0, b2=0, b3=0, alpha=0.1), trace=FALSE)

cat("fit fertility model m8g\n")
m8g <- mle2(y ~ dgeom(prob=1/(1+exp(b0 + b1*cen(dk$age) + b2*(cen(dk$age))^2 + b3*dk$phi_bar))), data=list(y = dk$n_kids), start=list(b0=0.2, b1=0, b2=0, b3=0), trace=FALSE)

cat("fit fertility model m9p\n")
m9p <- mle2(y ~ dpois(lambda=exp(b0 + b1*cen(dk$age) + b2*(cen(dk$age))^2 + b3*dk$snoob + b4*dk$phi_bar)), data=list(y = dk$n_kids), start=list(b0=0.2, b1=0.01, b2=0, b3=0, b4=0), trace=FALSE)

cat("fit fertility model m9z\n")
m9z <- mle2(y ~ dpois(lambda=exp(b0 + b1*cen(dk$age) + b2*(cen(dk$age))^2 + b3*dk$snoob + b4*dk$phi_bar)), data=list(y = dk$n_kids), start=list(b0=0.2, b1=0.01, b2=0, b3=0, b4=0), trace=FALSE)

cat("fit fertility model m9g\n")
m9g <- mle2(y ~ dpois(lambda=exp(b0 + b1*cen(dk$age) + b2*(cen(dk$age))^2 + b3*dk$snoob + b4*dk$phi_bar)), data=list(y = dk$n_kids), start=list(b0=0.2, b1=0.01, b2=0, b3=0, b4=0), trace=FALSE)

cat("fit fertility model m10p\n")
m10p <- mle2(y ~ dpois(lambda=exp(b0 + b1*cen(dk$age) + b2*(cen(dk$age))^2 + b3*dk$snoob + b4*dk$male + b5*dk$phi_bar)), data=list(y = dk$n_kids), start=list(b0=0.2, b1=0.01, b2=0, b3=0, b4=0, b5=0), trace=FALSE)

cat("fit fertility model m10z\n")
m10z <- mle2(y ~ my.dzpois(lambda=exp(b0 + b1*cen(dk$age) + b2*(cen(dk$age))^2 + b3*dk$snoob + b4*dk$male + b5*dk$phi_bar), a=alpha), data=list(y = dk$n_kids), start=list(b0=-1.15, b1=-0.02, b2=0, b3=0.68, b4=0.03, b5=0.22, alpha=0.1), trace=FALSE)

cat("fit fertility model m10g\n")
m10g <- mle2(y ~ dgeom(prob=1/(1+exp(b0 + b1*cen(dk$age) + b2*(cen(dk$age))^2 + b3*dk$snoob + b4*dk$male + b5*dk$phi_bar))), data=list(y = dk$n_kids), start=list(b0=0.2, b1=0.01, b2=0, b3=0, b4=0, b5=0), trace=FALSE)

AICtab(m10p, m10z, m10g, m9p, m9z, m9g, m8p, m8z, m8g, m7p, m7z, m7g, m4p, m4z, m4g, m3p, m3z, m3g, m2p, m2z, m2g, m1p, m1z, m1g, m0p, m0z, m0g, weights=T)

# interactions for the best models, m4 and m10

cat("fit fertility model m4p.i\n")
m4p.i <- mle2(y ~ dpois(lambda=exp(b0 + b1*cen(dk$age) + b2*(cen(dk$age))^2 + b3*dk$snoob + b4*dk$male + b5*dk$male*dk$snoob)), data=list(y = dk$n_kids), start=list(b0=0.2, b1=0.01, b2=0, b3=0, b4=0, b5=0), trace=FALSE)

cat("fit fertility model m4z.i\n")
m4z.i <- mle2(y ~ my.dzpois(lambda=exp(b0 + b1*cen(dk$age) + b2*(cen(dk$age))^2 + b3*dk$snoob + b4*dk$male + b5*dk$male*dk$snoob), a=alpha), data=list(y = dk$n_kids), start=list(b0=0.2, b1=0.01, b2=0, b3=0, b4=0, b5=0, alpha=0.1), trace=FALSE)

cat("fit fertility model m4g.i\n")
m4g.i <- mle2(y ~ dgeom(prob=1/(1+exp(b0 + b1*cen(dk$age) + b2*(cen(dk$age))^2 + b3*dk$snoob + b4*dk$male + b5*dk$male*dk$snoob))), data=list(y = dk$n_kids), start=list(b0=0.2, b1=0.01, b2=0, b3=0, b4=0, b5=0), trace=FALSE)

cat("fit fertility model m10p.i\n")
m10p.i <- mle2(y ~ dpois(lambda=exp(b0 + b1*cen(dk$age) + b2*(cen(dk$age))^2 + b3*dk$snoob + b4*dk$male + b5*dk$phi_bar + b6*dk$male*dk$snoob)), data=list(y = dk$n_kids), start=list(b0=0.2, b1=0.01, b2=0, b3=0, b4=0, b5=0, b6=0), trace=FALSE)

cat("fit fertility model m10z.i\n")
m10z.i <- mle2(y ~ my.dzpois(lambda=exp(b0 + b1*cen(dk$age) + b2*(cen(dk$age))^2 + b3*dk$snoob + b4*dk$male + b5*dk$phi_bar + b6*dk$male*dk$snoob), a=alpha), data=list(y = dk$n_kids), start=list(b0=-1.15, b1=-0.02, b2=0, b3=0.68, b4=0.03, b5=0.22, b6=0, alpha=0.1), trace=FALSE)

cat("fit fertility model m10g.i\n")
m10g.i <- mle2(y ~ dgeom(prob=1/(1+exp(b0 + b1*cen(dk$age) + b2*(cen(dk$age))^2 + b3*dk$snoob + b4*dk$male + b5*dk$phi_bar + b6*dk$male*dk$snoob))), data=list(y = dk$n_kids), start=list(b0=0.2, b1=0.01, b2=0, b3=0, b4=0, b5=0, b6=0), trace=FALSE)

AICtab(m10p.i, m10z.i, m10g.i, m4p.i, m4z.i, m4g.i, m10p, m10z, m10g, m9p, m9z, m9g, m8p, m8z, m8g, m7p, m7z, m7g, m4p, m4z, m4g, m3p, m3z, m3g, m2p, m2z, m2g, m1p, m1z, m1g, m0p, m0z, m0g, weights=T)

# this doesn't work for some reason
# compare(m10p.i, m10z.i, m10g.i, m4p.i, m4z.i, m4g.i, m10p, m10z, m10g, m9p, m9z, m9g, m8p, m8z, m8g, m7p, m7z, m7g, m4p, m4z, m4g, m3p, m3z, m3g, m2p, m2z, m2g, m1p, m1z, m1g, m0p, m0z, m0g, nobs=length(dk$n_kids))

# this figure requires snoobs to be near 50% prevalence at some point
if (any(abs(cens$phi_bar - 0.5) < 0.1)) {

  b0 <- -0.878
  b1 <- -0.006
  b2 <- -0.002
  b3 <- 0.681
  b4 <- 0.120
  b5 <- 0.244
  alpha <- 0.357

  cat("create figure fertility.png\n")

  png("./figures/fertility.png", height=15, width=20, 
          units='cm', res=300)

  par(mfrow=c(2,2))

  age_cen <- 32.67694

  k <- 1
  plot.phibar <- 0.5
  will.have.k <- as.numeric(dk$n_kids==k)
  dk$age[dk$age >= 60] <- 59
  s.pr.will.have.k <- tapply(will.have.k[dk$snoob==1 & abs(dk$phi_bar - plot.phibar) < 0.1], floor(dk$age[dk$snoob==1 & abs(dk$phi_bar - 0.5) < 0.1]), mean)
  ns.pr.will.have.k <- tapply(will.have.k[dk$snoob==0 & abs(dk$phi_bar - plot.phibar) < 0.1], floor(dk$age[dk$snoob==0 & abs(dk$phi_bar - 0.5) < 0.1]), mean)
  age.list <- 11:59

  plot(ns.pr.will.have.k ~ age.list, ylim=c(min(c(ns.pr.will.have.k, s.pr.will.have.k)), max(ns.pr.will.have.k, s.pr.will.have.k)), ylab="", xlab="Age", las=1, main=paste("Prob. of having", k, "kid", sep=" "), col="gray")
  points(s.pr.will.have.k ~ age.list, col="black", pch=1)

  curve(my.dzpois(k, lambda=exp(b0 + b1*(x-age_cen) + b2*(x-age_cen)^2 + b3*1 + b4*1 + b5*plot.phibar), a=alpha), from=11, to=60, col="dodgerblue", ylim=c(0, 0.2), add=T)
  curve(my.dzpois(k, lambda=exp(b0 + b1*(x-age_cen) + b2*(x-age_cen)^2 + b3*0 + b4*1 + b5*plot.phibar), a=alpha), from=11, to=60, col="dodgerblue", lty=1, add=T)
  curve(my.dzpois(k, lambda=exp(b0 + b1*(x-age_cen) + b2*(x-age_cen)^2 + b3*1 + b4*0 + b5*plot.phibar), a=alpha), from=11, to=60, col="maroon1", add=T)
  curve(my.dzpois(k, lambda=exp(b0 + b1*(x-age_cen) + b2*(x-age_cen)^2 + b3*0 + b4*0 + b5*plot.phibar), a=alpha), from=11, to=60, col="maroon1", lty=1, add=T)

  k <- 2
  plot.phibar <- 0.5
  will.have.k <- as.numeric(dk$n_kids==k)
  s.pr.will.have.k <- tapply(will.have.k[dk$snoob==1 & abs(dk$phi_bar - plot.phibar) < 0.1], floor(dk$age[dk$snoob==1 & abs(dk$phi_bar - 0.5) < 0.1]), mean)
  ns.pr.will.have.k <- tapply(will.have.k[dk$snoob==0 & abs(dk$phi_bar - plot.phibar) < 0.1], floor(dk$age[dk$snoob==0 & abs(dk$phi_bar - 0.5) < 0.1]), mean)
  age.list <- 11:59

  plot(ns.pr.will.have.k ~ age.list, ylim=c(min(c(ns.pr.will.have.k, s.pr.will.have.k)), max(ns.pr.will.have.k, s.pr.will.have.k)), ylab="", xlab="Age", las=1, main=paste("Prob. of having", k, "kids", sep=" "), col="gray")
  points(s.pr.will.have.k ~ age.list, col="black", pch=1)

  curve(my.dzpois(k, lambda=exp(b0 + b1*(x-age_cen) + b2*(x-age_cen)^2 + b3*1 + b4*1 + b5*plot.phibar), a=alpha), from=11, to=60, col="dodgerblue", ylim=c(0, 0.2), add=T)
  curve(my.dzpois(k, lambda=exp(b0 + b1*(x-age_cen) + b2*(x-age_cen)^2 + b3*0 + b4*1 + b5*plot.phibar), a=alpha), from=11, to=60, col="dodgerblue", lty=1, add=T)
  curve(my.dzpois(k, lambda=exp(b0 + b1*(x-age_cen) + b2*(x-age_cen)^2 + b3*1 + b4*0 + b5*plot.phibar), a=alpha), from=11, to=60, col="maroon1", add=T)
  curve(my.dzpois(k, lambda=exp(b0 + b1*(x-age_cen) + b2*(x-age_cen)^2 + b3*0 + b4*0 + b5*plot.phibar), a=alpha), from=11, to=60, col="maroon1", lty=1, add=T)


  k <- 3
  plot.phibar <- 0.5
  will.have.k <- as.numeric(dk$n_kids==k)
  s.pr.will.have.k <- tapply(will.have.k[dk$snoob==1 & abs(dk$phi_bar - plot.phibar) < 0.1], floor(dk$age[dk$snoob==1 & abs(dk$phi_bar - 0.5) < 0.1]), mean)
  ns.pr.will.have.k <- tapply(will.have.k[dk$snoob==0 & abs(dk$phi_bar - plot.phibar) < 0.1], floor(dk$age[dk$snoob==0 & abs(dk$phi_bar - 0.5) < 0.1]), mean)
  age.list <- 11:59

  plot(ns.pr.will.have.k ~ age.list, ylim=c(min(c(ns.pr.will.have.k, s.pr.will.have.k)), max(ns.pr.will.have.k, s.pr.will.have.k)), ylab="", xlab="Age", las=1, main=paste("Prob. of having", k, "kids", sep=" "), col="gray")
  points(s.pr.will.have.k ~ age.list, col="black", pch=1)

  curve(my.dzpois(k, lambda=exp(b0 + b1*(x-age_cen) + b2*(x-age_cen)^2 + b3*1 + b4*1 + b5*plot.phibar), a=alpha), from=11, to=60, col="dodgerblue", ylim=c(0, 0.2), add=T)
  curve(my.dzpois(k, lambda=exp(b0 + b1*(x-age_cen) + b2*(x-age_cen)^2 + b3*0 + b4*1 + b5*plot.phibar), a=alpha), from=11, to=60, col="dodgerblue", lty=1, add=T)
  curve(my.dzpois(k, lambda=exp(b0 + b1*(x-age_cen) + b2*(x-age_cen)^2 + b3*1 + b4*0 + b5*plot.phibar), a=alpha), from=11, to=60, col="maroon1", add=T)
  curve(my.dzpois(k, lambda=exp(b0 + b1*(x-age_cen) + b2*(x-age_cen)^2 + b3*0 + b4*0 + b5*plot.phibar), a=alpha), from=11, to=60, col="maroon1", lty=1, add=T)


  k <- 4
  plot.phibar <- 0.5
  will.have.k <- as.numeric(dk$n_kids==k)
  s.pr.will.have.k <- tapply(will.have.k[dk$snoob==1 & abs(dk$phi_bar - plot.phibar) < 0.1], floor(dk$age[dk$snoob==1 & abs(dk$phi_bar - 0.5) < 0.1]), mean)
  ns.pr.will.have.k <- tapply(will.have.k[dk$snoob==0 & abs(dk$phi_bar - plot.phibar) < 0.1], floor(dk$age[dk$snoob==0 & abs(dk$phi_bar - 0.5) < 0.1]), mean)
  age.list <- 11:59

  plot(ns.pr.will.have.k ~ age.list, ylim=c(min(c(ns.pr.will.have.k, s.pr.will.have.k)), max(ns.pr.will.have.k, s.pr.will.have.k)), ylab="", xlab="Age", las=1, main=paste("Prob. of having", k, "kids", sep=" "), col="gray")
  points(s.pr.will.have.k ~ age.list, col="black", pch=1)

  curve(my.dzpois(k, lambda=exp(b0 + b1*(x-age_cen) + b2*(x-age_cen)^2 + b3*1 + b4*1 + b5*plot.phibar), a=alpha), from=11, to=60, col="dodgerblue", ylim=c(0, 0.2), add=T)
  curve(my.dzpois(k, lambda=exp(b0 + b1*(x-age_cen) + b2*(x-age_cen)^2 + b3*0 + b4*1 + b5*plot.phibar), a=alpha), from=11, to=60, col="dodgerblue", lty=1, add=T)
  curve(my.dzpois(k, lambda=exp(b0 + b1*(x-age_cen) + b2*(x-age_cen)^2 + b3*1 + b4*0 + b5*plot.phibar), a=alpha), from=11, to=60, col="maroon1", add=T)
  curve(my.dzpois(k, lambda=exp(b0 + b1*(x-age_cen) + b2*(x-age_cen)^2 + b3*0 + b4*0 + b5*plot.phibar), a=alpha), from=11, to=60, col="maroon1", lty=1, add=T)

  dev.off()

}



# Individual Change
cat("fit individual change models\n")

cat("fit change model m0\n")

m0 <- mle2(y ~ dbinom(1, prob = p), data = list(y = d$snoob), start = list(p = 0.2), trace = FALSE)

cat("fit change model m1\n")

m1 <- mle2(y ~ dbinom(1, prob=L*(x^(B)/(x^(B) + (1-x)^(B))) + (1-L)*z),
  data = list(y = d$will_be_snoob, x = d$phi_bar, z = d$snoob), start = list(L = 0.2, B = 0.2), trace = FALSE) # warnings produced

cat("fit change model m2\n")

m2 <- mle2(y ~ dbinom(1, prob=L*(x + 4*(B-1)*x*(1-x)*(x-0.5))+ (1-L)*z), data=list(y = d$will_be_snoob, x =d$phi_bar, z=d$snoob), start=list(L=0.2, B=0), trace=FALSE)

cat("fit change model m3\n")

m3 <- mle2(y ~ dbinom(1, prob=L*(x + 4*(B-1)*x*(1-x)*(x-k))+ (1-L)*z), data=list(y = d$will_be_snoob, x =d$phi_bar, z=d$snoob), start=list(L=0.2, B=0, k=0.5), trace=FALSE)

cat("fit change model m4\n")

m4 <-  mle2(y ~ dbinom(1, prob=L*(exp(b0 + b1*x)/(1+exp(b0 + b1*x))) + (1-L)*z ), data=list(y = d$will_be_snoob, x = d$phi_bar, z = d$snoob), start=list(L=0.2, b0=0.3, b1=0), trace=FALSE)

cat("fit change model m5\n")

m5 <- mle2(y ~ dbinom(1, prob=L*(1/(1+exp(b0 + b1*n))) + (1-L)*z), data=list(y = d$will_be_snoob, n = d$pop_size, z = d$snoob), start=list(L=0.2, b0=0.3, b1=0), trace=FALSE)

cat("fit change model m6\n")

m6 <- mle2(y ~ dbinom(1, prob=L*(1/(1+exp(b0 + b1*su(n) + b2*su(n)^2))) + (1-L)*z), data=list(y = d$will_be_snoob, n = d$pop_size, z = d$snoob), start=list(L=0.2, b0=0.3, b1=0, b2=0), trace=FALSE)

cat("fit change model m7\n")

m7 <- mle2(y ~ dbinom(1, prob=L*(1/(1+exp(b0 + b1*t))) + (1-L)*z ), data=list(y = d$will_be_snoob, t = d$census, z = d$snoob), start=list(L=0.2, b0=0.3, b1=0), trace=FALSE)

cat("fit change model m8\n")

m8 <- mle2(y ~ dbinom(1, prob=L*(1/(1+exp(b0 + b1*t + b2*t^2))) + (1-L)*z ), data=list(y = d$will_be_snoob, t = d$census, z = d$snoob), start=list(L=0.2, b0=0.3, b1=0, b2=0), trace=FALSE)

cat("fit change model m9\n")

m9 <- mle2(y ~ dbinom(1, prob= L*(b0 + b1*sin(b2*t)) + (1-L)*z), data=list(y = d$will_be_snoob, t = d$census, z = d$snoob), start=list(L=0.2, b0=0.3, b1=0.1, b2=1), trace=FALSE)

cat("fit change model m10\n")

m10 <- mle2(y ~ dbinom(1, prob= L*(b0 + b1*sin(b2*n)) + (1-L)*z), data=list(y = d$will_be_snoob, n = d$pop_size, z = d$snoob), start=list(L=0.2, b0=0.3, b1=0.1, b2=1), trace=FALSE)

cat("fit change model m11\n")

m11 <- mle2(y ~ dbinom(1, prob=L*(1/(1+exp(b0 + b1*n + b2*t))) + (1-L)*z), data=list(y = d$will_be_snoob, n = d$pop_size, t = d$census, z = d$snoob), start=list(L=0.2, b0=0.3, b1=0, b2=0), trace=FALSE)

cat("fit change model m12\n")

m12 <-  mle2(y ~ dbinom(1, prob=L*(1/(1+exp(b0 + b1*x + b2*a + b3*m))) + (1-L)*z ), data=list(y = d$will_be_snoob, x = d$phi_bar, z = d$snoob, a = cen(d$age), m = d$male), start=list(L=0.2, b0=0.3, b1=0, b2=0, b3=0), trace=FALSE)

cat("fit change model m13\n")

m13 <-  mle2(y ~ dbinom(1, prob=L*(1/(1+exp(b0 + b1*x + b2*a + b3*m + b4*t))) + (1-L)*z ), data=list(y = d$will_be_snoob, x = d$phi_bar, a = cen(d$age), m = d$male, z = d$snoob, t = d$census), start=list(L=0.2, b0=0.3, b1=0, b2=0, b3=0, b4=0), trace=FALSE)

cat("create figure learning.png\n")

png("./figures/learning.png", height=15, width=20, units='cm', res=300)

cfact <- data.frame(
  frq_snoob = seq(0.01, 0.99, length=100),
  pr_stay_snoob_mu = NA,
  pr_stay_snoob_lb = NA,
  pr_stay_snoob_ub = NA,
  pr_to_snoob_mu = NA,
  pr_to_snoob_lb = NA,
  pr_to_snoob_ub = NA
)

post <- rethinking::extract.samples(m1)

for (i in seq_len(nrow(cfact))) {
  prob <- post$L * (cfact$frq_snoob[i])^(post$B)/(cfact$frq_snoob[i]^post$B + (1 - cfact$frq_snoob[i])^post$B)
  cfact$pr_to_snoob_mu[i] <- mean(prob)
  cfact$pr_to_snoob_lb[i] <- HPDI(prob)[1]
  cfact$pr_to_snoob_ub[i] <- HPDI(prob)[2]
  prob <- post$L * (cfact$frq_snoob[i])^(post$B)/(cfact$frq_snoob[i]^post$B + (1 - cfact$frq_snoob[i])^post$B) + (1 - post$L)
  cfact$pr_stay_snoob_mu[i] <- mean(prob)
  cfact$pr_stay_snoob_lb[i] <- HPDI(prob)[1]
  cfact$pr_stay_snoob_ub[i] <- HPDI(prob)[2]
}

par(mfrow = c(1, 2))

plot(cens$snoob_to_snoob ~ cens$phi_bar,
  xlim = c(0.15, 1), ylim = c(0.85, 1),
  ylab = "Prob. Remaining Snoob", xlab="Mean Pop. Frq",
  col = "blue", las = 1, frame.plot = FALSE)

points(cfact$frq_snoob, cfact$pr_stay_snoob_mu, lty = 2, type = "l")
points(cfact$frq_snoob, cfact$pr_stay_snoob_lb, lty = 1, type = "l")
points(cfact$frq_snoob, cfact$pr_stay_snoob_ub, lty = 1, type = "l")

plot(cens$non_to_snoob ~ cens$phi_bar,
  xlim = c(0.15, 1), ylim = c(0, 0.15),
  ylab = "Prob. Becoming Snoob", xlab = "Mean Pop. Frq",
  col = "blue", las = 1, frame.plot = FALSE)

points(cfact$frq_snoob, cfact$pr_to_snoob_mu, lty = 2, type = "l")
points(cfact$frq_snoob, cfact$pr_to_snoob_lb, lty = 1, type = "l")
points(cfact$frq_snoob, cfact$pr_to_snoob_ub, lty = 1, type = "l")

dev.off()
