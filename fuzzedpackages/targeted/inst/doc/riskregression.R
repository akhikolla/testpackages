## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
library(lava)

## ---- fig.width=5, fig.height=5-----------------------------------------------
  p0 <- seq(0,1,length.out=100)
  p1 <- function(p0,op) 1/(1+(op*(1-p0)/p0)^-1)
  plot(0, type="n", xlim=c(0,1), ylim=c(0,1),
     xlab="P(Y=1|A=0)", ylab="P(Y=1|A=1)", main="Constant odds product")
  for (op in exp(seq(-6,6,by=.25))) lines(p0,p1(p0,op), col="lightblue")

## ---- fig.width=5, fig.height=5, fig.caption=' '------------------------------
  p0 <- seq(0,1,length.out=100)
  p1 <- function(p0,rr) rr*p0
  plot(0, type="n", xlim=c(0,1), ylim=c(0,1),
     xlab="P(Y=1|A=0)", ylab="P(Y=1|A=1)", main="Constant relative risk")
  for (rr in exp(seq(-3,3,by=.25))) lines(p0,p1(p0,rr), col="lightblue")

## ----setup--------------------------------------------------------------------
library(targeted)

## -----------------------------------------------------------------------------
m <- lvm(a ~ x,
         lp.target ~ 1,
         lp.nuisance ~ x+z)
m <- binomial.rr(m, response="y", exposure="a", target.model="lp.target", nuisance.model="lp.nuisance")

## -----------------------------------------------------------------------------
args(binomial.rr)

## -----------------------------------------------------------------------------
args(binomial.rd)

## -----------------------------------------------------------------------------
coef(m)

## -----------------------------------------------------------------------------
p <- c('a'=-1, 'lp.target'=1, 'lp.nuisance'=-1, 'a~x'=2)

## -----------------------------------------------------------------------------
d <- sim(m, 1e4, p=p, seed=1)

head(d)

## -----------------------------------------------------------------------------
args(riskreg_mle)

## -----------------------------------------------------------------------------
x1 <- model.matrix(~1, d)
x2 <- model.matrix(~x+z, d)

fit1 <- with(d, riskreg_mle(y, a, x1, x2, type="rr"))
fit1

## -----------------------------------------------------------------------------
estimate(fit1, keep=1)

## -----------------------------------------------------------------------------
with(d, riskreg_fit(y, a, target=x1, nuisance=x2, propensity=x2, type="rr"))


## -----------------------------------------------------------------------------
  riskreg(y~a, nuisance=~x+z, propensity=~z, data=d, type="rr")

## -----------------------------------------------------------------------------
  riskreg(y~a, nuisance=~z, propensity=~x+z, data=d, type="rr")

## -----------------------------------------------------------------------------
  fit2 <- with(d, riskreg_mle(y, a, x1=model.matrix(~1,d), x2=model.matrix(~z, d)))
  estimate(fit2, keep=1)

## -----------------------------------------------------------------------------
fit <- riskreg(y~a, target=~x, nuisance=~x+z, data=d)
fit

## -----------------------------------------------------------------------------
m2 <- binomial.rd(m, response="y", exposure="a", target.model="lp.target", nuisance.model="lp.nuisance")
d2 <- sim(m2, 1e4, p=p)

## -----------------------------------------------------------------------------
riskreg(y~a, nuisance=~x+z, data=d2, type="rd")

## -----------------------------------------------------------------------------
head(iid(fit))

## -----------------------------------------------------------------------------
sessionInfo()

