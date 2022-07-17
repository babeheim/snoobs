
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
