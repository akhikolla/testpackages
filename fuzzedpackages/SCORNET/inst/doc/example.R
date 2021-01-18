## ----setup---------------------------------------------------------------
library(SCORNET)

## ------------------------------------------------------------------------
sim <- function(N){
  dat <- data.frame('ID'=1:N)
  dat$Z0 <- runif(N,-1,1)
  dat$C <- 10*(rexp(N)*exp(-0.75*dat$Z))^(2/3)
  dat$T <- 15*(rexp(N)*exp(-0.75*dat$Z))^(2/5)
  dat$X <- pmin(dat$T,dat$C)
  dat$Delta <- dat$T <= dat$C
  dat$filter <- as.logical(rbinom(N,1,0.98)*dat$Delta + rbinom(N,1,0.12)*(1-dat$Delta))
  dat$Zehr <- pmin(dat$T+rnorm(N,0,2),dat$C)
  return(dat)
}

N <- 10000
n <- 200
dat <- sim(N)

## ---- message=FALSE------------------------------------------------------
t0.all <- seq(quantile(dat$C,.1),quantile(dat$C,.9),length.out=100)

t1 <- proc.time()
scornet_est <- scornet(dat$Delta[1:n],dat$C[1:n],t0.all,dat$C[-c(1:n)],dat$filter[1:n],dat$filter[-c(1:n)],dat$Z0[1:n],dat$Z0[-c(1:n)],dat$Zehr[1:n],dat$Zehr[-c(1:n)])
t2 <- proc.time()

## ---- message=FALSE------------------------------------------------------
plot(t0.all,sapply(t0.all,function(t){mean(dat$T>t)}),xlab='Time',ylab='Survival',
col='black',type='l')
lines(t0.all,scornet_est$S_hat,col='red',type='l')

## ---- message=FALSE------------------------------------------------------
plot(t0.all,scornet_est$StdErrs,type='l',xlab='Time',ylab='Wald Std. Err.')

## ---- message=FALSE------------------------------------------------------
print(t2-t1)

