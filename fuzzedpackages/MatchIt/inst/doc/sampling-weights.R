## ---- include = FALSE---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
knitr::opts_chunk$set(echo = TRUE, eval=T)
options(width = 200, digits = 4)

#Generatng data similar to Austin (2009) for demonstrating treatment effect estimation with sampling weights
gen_X <- function(n) {
  X <- matrix(rnorm(9 * n), nrow = n, ncol = 9)
  X[,5] <- as.numeric(X[,5] < .5)
  X
}

#~20% treated
gen_A <- function(X) {
  LP_A <- - 1.2 + log(2)*X[,1] - log(1.5)*X[,2] + log(2)*X[,4] - log(2.4)*X[,5] + log(2)*X[,7] - log(1.5)*X[,8]
  P_A <- plogis(LP_A)
  rbinom(nrow(X), 1, P_A)
}

# Continuous outcome
gen_Y_C <- function(A, X) {
  2*A + 2*X[,1] + 2*X[,2] + 2*X[,3] + 1*X[,4] + 2*X[,5] + 1*X[,6] + rnorm(length(A), 0, 5)
}
#Conditional:
#  MD: 2
#Marginal:
#  MD: 2

gen_SW <- function(X) {
  e <- rbinom(nrow(X), 1, .3)
  1/plogis(log(1.4)*X[,2] + log(.7)*X[,4] + log(.9)*X[,6] + log(1.5)*X[,8] + log(.9)*e +
             -log(.5)*e*X[,2] + log(.6)*e*X[,4])
}

set.seed(19599)

n <- 2000
X <- gen_X(n)
A <- gen_A(X)
SW <- gen_SW(X)

Y_C <- gen_Y_C(A, X)

d <- data.frame(A, X, Y_C, SW)

## ----message=FALSE,warning=FALSE----------------------------------------------------------------------------------------------------------------------------------------------------------------------
head(d)

library("MatchIt")
library("lmtest")
library("sandwich")

## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
mF_s <- matchit(A ~ X1 + X2 + X3 + X4 + X5 + 
                  X6 + X7 + X8 + X9, data = d,
                method = "full", distance = "glm",
                estimand = "ATE", s.weights = ~SW)
mF_s

## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
mF <- matchit(A ~ X1 + X2 + X3 + X4 + X5 + 
                X6 + X7 + X8 + X9, data = d,
              method = "full", distance = "glm",
              estimand = "ATE")
mF

## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
mF <- add_s.weights(mF, ~SW)

mF

## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Balance before matching and for the SW propensity score full matching
summary(mF_s, improvement = FALSE)

#Balance for the non-SW propensity score full matching
summary(mF, un = FALSE)

## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
library(sandwich)
library(lmtest)

md_F_s <- match.data(mF_s)

fit <- lm(Y_C ~ A + X1 + X2 + X3 + X4 + X5 + 
             X6 + X7 + X8 + X9, data = md_F_s,
          weights = weights)

coeftest(fit, vcov. = vcovCL, cluster = ~subclass)["A",,drop=FALSE]

## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
md_F_s <- match.data(mF_s, include.s.weights = FALSE)

fit <- lm(Y_C ~ A + X1 + X2 + X3 + X4 + X5 + 
             X6 + X7 + X8 + X9, data = md_F_s,
          weights = weights * SW)

coeftest(fit, vcov. = vcovCL, cluster = ~subclass)["A",,drop=FALSE]

## ---- eval = FALSE------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#  #Generatng data similar to Austin (2009) for demonstrating
#  #treatment effect estimation with sampling weights
#  gen_X <- function(n) {
#    X <- matrix(rnorm(9 * n), nrow = n, ncol = 9)
#    X[,5] <- as.numeric(X[,5] < .5)
#    X
#  }
#  
#  #~20% treated
#  gen_A <- function(X) {
#    LP_A <- - 1.2 + log(2)*X[,1] - log(1.5)*X[,2] + log(2)*X[,4] - log(2.4)*X[,5] +
#      log(2)*X[,7] - log(1.5)*X[,8]
#    P_A <- plogis(LP_A)
#    rbinom(nrow(X), 1, P_A)
#  }
#  
#  # Continuous outcome
#  gen_Y_C <- function(A, X) {
#    2*A + 2*X[,1] + 2*X[,2] + 2*X[,3] + 1*X[,4] + 2*X[,5] + 1*X[,6] + rnorm(length(A), 0, 5)
#  }
#  #Conditional:
#  #  MD: 2
#  #Marginal:
#  #  MD: 2
#  
#  gen_SW <- function(X) {
#    e <- rbinom(nrow(X), 1, .3)
#    1/plogis(log(1.4)*X[,2] + log(.7)*X[,4] + log(.9)*X[,6] + log(1.5)*X[,8] + log(.9)*e +
#               -log(.5)*e*X[,2] + log(.6)*e*X[,4])
#  }
#  
#  set.seed(19599)
#  
#  n <- 2000
#  X <- gen_X(n)
#  A <- gen_A(X)
#  SW <- gen_SW(X)
#  
#  Y_C <- gen_Y_C(A, X)
#  
#  d <- data.frame(A, X, Y_C, SW)

