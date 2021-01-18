## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)
## Added 20190717 to get vignette to build
pkgbuild::compile_dll()

## ---- candlestick-------------------------------------------------------------
# candlestick function
# J C Nash 2011-2-3
cstick.f<-function(x,alpha=100){
  x<-as.vector(x)
  r2<-crossprod(x)
  f<-as.double(r2+alpha/r2)
  return(f)
}

cstick.g<-function(x,alpha=100){
  x<-as.vector(x)
  r2<-as.numeric(crossprod(x))
  g1<-2*x
  g2 <- (-alpha)*2*x/(r2*r2)
  g<-as.double(g1+g2)
  return(g)
}
library(lbfgsb3c)
nn <- 2
x0 <- c(10,10)
lo <- c(1, 1)
up <- c(10,10)
print(x0)
## c2o <- opm(x0, cstick.f, cstick.g, lower=lo, upper=up, method=meths, control=list(trace=0))
## print(summary(c2o, order=value))
c2l1 <- lbfgsb3c(x0, cstick.f, cstick.g, lower=lo, upper=up)
c2l1

## meths <- c("L-BFGS-B", "lbfgsb3c", "Rvmmin", "Rcgmin", "Rtnmin")
## require(optimx)

## cstick2a <- opm(x0, cstick.f, cstick.g, method=meths, upper=up, lower=lo, control=list(kkt=FALSE))
## print(summary(cstick2a, par.select=1:2, order=value))
lo <- c(4, 4)
## c2ob <- opm(x0, cstick.f, cstick.g, lower=lo, upper=up, method=meths, control=list(trace=0))
## print(summary(c2ob, order=value))
c2l2 <- lbfgsb3c(x0, cstick.f, cstick.g, lower=lo, upper=up)
c2l2

## cstick2b <- opm(x0, cstick.f, cstick.g, method=meths, upper=up, lower=lo, control=list(kkt=FALSE))
## print(summary(cstick2b, par.select=1:2, order=value))

## nn <- 100
## x0 <- rep(10, nn)
## up <- rep(10, nn)
## lo <- rep(1e-4, nn)
## cco <- opm(x0, cstick.f, cstick.g, lower=lo, upper=up, method=meths, control=list(trace=0, kkt=FALSE))
## print(summary(cco, par.select=1:4, order=value))

## ----exrosen------------------------------------------------------------------
# require(funconstrain) ## not in CRAN, so explicit inclusion of this function
# exrosen <- ex_rosen()
# exrosenf <- exrosen$fn
exrosenf <- function (par) {
    n <- length(par)
    if (n%%2 != 0) {
        stop("Extended Rosenbrock: n must be even")
    }
    fsum <- 0
    for (i in 1:(n/2)) {
        p2 <- 2 * i
        p1 <- p2 - 1
        f_p1 <- 10 * (par[p2] - par[p1]^2)
        f_p2 <- 1 - par[p1]
        fsum <- fsum + f_p1 * f_p1 + f_p2 * f_p2
    }
    fsum
}
# exroseng <- exrosen$gr
exroseng <- function (par) {
    n <- length(par)
    if (n%%2 != 0) {
        stop("Extended Rosenbrock: n must be even")
    }
    grad <- rep(0, n)
    for (i in 1:(n/2)) {
        p2 <- 2 * i
        p1 <- p2 - 1
        xx <- par[p1] * par[p1]
        yx <- par[p2] - xx
        f_p1 <- 10 * yx
        f_p2 <- 1 - par[p1]
        grad[p1] <- grad[p1] - 400 * par[p1] * yx - 2 * f_p2
        grad[p2] <- grad[p2] + 200 * yx
    }
    grad
}

exrosenx0 <- function (n = 20) {
    if (n%%2 != 0) {
        stop("Extended Rosenbrock: n must be even")
    }
    rep(c(-1.2, 1), n/2)
}


require(lbfgsb3c)
## require(optimx)

## require(optimx)
for (n in seq(2,12, by=2)) {
  cat("ex_rosen try for n=",n,"\n")
  x0 <- exrosenx0(n)
  lo <- rep(-1.5, n)
  up <- rep(3, n)
  print(x0)
  cat("optim L-BFGS-B\n")
  eo <- optim(x0, exrosenf, exroseng, lower=lo, upper=up, method="L-BFGS-B", control=list(trace=0))
  print(eo)
  cat("lbfgsb3c\n")
  el <- lbfgsb3c(x0, exrosenf, exroseng, lower=lo, upper=up, control=list(trace=0))
  print(el)
##    erfg <- opm(x0, exrosenf, exroseng, method=meths, lower=lo, upper=up)
##    print(summary(erfg, par.select=1:2, order=value))
}

## ---- usingFortran, eval=FALSE------------------------------------------------
#  system("R CMD SHLIB jrosen.f")
#  dyn.load("jrosen.so")
#  is.loaded("rosen")
#  x0 <- as.double(c(-1.2,1))
#  fv <- as.double(-999)
#  n <- as.double(2)
#  testf <- .Fortran("rosen", n=as.integer(n), x=as.double(x0), fval=as.double(fv))
#  testf
#  
#  rrosen <- function(x) {
#    fval <- 0.0
#    for (i in 1:(n-1)) {
#      dx <- x[i + 1] - x[i] * x[i]
#      fval <- fval + 100.0 * dx * dx
#      dx <- 1.0 - x[i]
#      fval <- fval + dx * dx
#    }
#    fval
#  }
#  
#  (rrosen(x0))
#  
#  frosen <- function(x){
#    nn <- length(x)
#    if (nn > 100) { stop("max number of parameters is 100")}
#    fv <- -999.0
#    val <- .Fortran("rosen", n=as.integer(nn), x=as.double(x), fval=as.double(fv))
#    val$fval # NOTE--need ONLY function value returned
#  }
#  # Test the funcion
#  tval <- frosen(x0)
#  str(tval)
#  
#  cat("Run with Nelder-Mead using R function\n")
#  mynm <- optim(x0, rrosen, control=list(trace=0))
#  print(mynm)
#  cat("\n\n Run with Nelder-Mead using Fortran function")
#  mynmf <- optim(x0, frosen, control=list(trace=0))
#  print(mynmf)
#  
#  
#  library(lbfgsb3c)
#  library(microbenchmark)
#  cat("try lbfgsb3c, no Gradient \n")
#  cat("R function\n")
#  tlR<-microbenchmark(myopR <- lbfgsb3c(x0, rrosen, gr=NULL, control=list(trace=0)))
#  print(tlR)
#  print(myopR)
#  cat("Fortran function\n")
#  tlF<-microbenchmark(myop <- lbfgsb3c(x0, frosen, gr=NULL, control=list(trace=0)))
#  print(tlF)
#  print(myop)

