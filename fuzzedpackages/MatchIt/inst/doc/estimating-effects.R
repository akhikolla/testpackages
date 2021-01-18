## ---- include = FALSE---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
knitr::opts_chunk$set(echo = TRUE, eval=T)
options(width = 200, digits= 4)

#Generatng data similar to Austin (2009) for demonstrating treatment effect estimation
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

# Binary outcome
gen_Y_B <- function(A, X) {
  LP_B <- -2 + log(2.4)*A + log(2)*X[,1] + log(2)*X[,2] + log(2)*X[,3] + log(1.5)*X[,4] + log(2.4)*X[,5] + log(1.5)*X[,6]
  P_B <- plogis(LP_B)
  rbinom(length(A), 1, P_B)
}
#Conditional:
#  OR:   2.4
#  logOR: .875
#Marginal:
#  RD:    .144
#  RR:   1.54
#  logRR: .433
#  OR:   1.92
#  logOR  .655

# Survival outcome
gen_Y_S <- function(A, X) {
  LP_S <- -2 + log(2.4)*A + log(2)*X[,1] + log(2)*X[,2] + log(2)*X[,3] + log(1.5)*X[,4] + log(2.4)*X[,5] + log(1.5)*X[,6]
  sqrt(-log(runif(length(A)))*2e4*exp(-LP_S))
}
#Conditional:
#  HR:   2.4
#  logHR: .875
#Marginal:
#  HR:   1.57
#  logHR: .452

set.seed(19599)

n <- 2000
X <- gen_X(n)
A <- gen_A(X)

Y_C <- gen_Y_C(A, X)
Y_B <- gen_Y_B(A, X)
Y_S <- gen_Y_S(A, X)

d <- data.frame(A, X, Y_C, Y_B, Y_S)

## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
head(d)

## ----message=FALSE,warning=FALSE----------------------------------------------------------------------------------------------------------------------------------------------------------------------
library("MatchIt")
library("lmtest")
library("sandwich")
library("boot")
library("survival")

## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
mNN <- matchit(A ~ X1 + X2 + X3 + X4 + X5 + 
                 X6 + X7 + X8 + X9, data = d)

mNN

md <- match.data(mNN)

head(md)

## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Linear model without covariates
fit1 <- lm(Y_C ~ A, data = md, weights = weights)

## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Cluster-robust standard errors
coeftest(fit1, vcov. = vcovCL, cluster = ~subclass)

## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Block bootstrap confidence interval
# library(boot)

pair_ids <- levels(md$subclass)

est_fun <- function(pairs, i) {
  
  #Compute number of times each pair is present
  numreps <- table(pairs[i])
  
  #For each pair p, copy corresponding md row indices numreps[p] times
  ids <- unlist(lapply(pair_ids[pair_ids %in% names(numreps)],
                       function(p) rep(which(md$subclass == p), 
                                              numreps[p])))
  
  #Subset md with block bootstrapped ids
  md_boot <- md[ids,]
  
  #Effect estimation
  fit_boot <- lm(Y_C ~ A, data = md_boot, weights = weights)
  
  #Return the coefficient on treatment
  return(coef(fit_boot)["A"])
}

boot_est <- boot(pair_ids, est_fun, R = 499)
boot_est
boot.ci(boot_est, type = "bca")

## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Linear model with covariates
fit3 <- lm(Y_C ~ A + X1 + X2 + X3 + X4 + X5 + 
             X6 + X7 + X8 + X9, data = md,
           weights = weights)

coeftest(fit3, vcov. = vcovCL, cluster = ~subclass)["A",,drop=FALSE]

## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Generalized linear model without covariates
fit4 <- glm(Y_B ~ A, data = md, family = binomial(link = "logit"))

## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Cluster-robust standard errors
coeftest(fit4, vcov. = vcovCL, cluster = ~subclass)
exp(coef(fit4)) #OR

## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Block bootstrap confidence interval
# library(boot)

pair_ids <- levels(md$subclass)

est_fun <- function(pairs, i) {
  
  #Compute number of times each pair is present
  numreps <- table(pairs[i])
  
  #For each pair p, copy corresponding md row indices numreps[p] times
  ids <- unlist(lapply(pair_ids[pair_ids %in% names(numreps)],
                       function(p) rep(which(md$subclass == p), 
                                              numreps[p])))
  
  #Subset md with block bootstrapped ids
  md_boot <- md[ids,]
  
  #Fitting outcome the model
  fit_boot <- glm(Y_B ~ A + X1 + X2 + X3 + X4 + X5 + 
                    X6 + X7 + X8 + X9, data = md_boot, 
                  family = binomial(link = "logit"),
                  weights = weights)
  
  #Estimate potential outcomes for each unit
  #Under control
  md_boot$A <- 0
  P0 <- weighted.mean(predict(fit_boot, md_boot, type = "response"),
                      w = md_boot$weights)
  Odds0 <- P0 / (1 - P0)
  
  #Under treatment
  md_boot$A <- 1
  P1 <- weighted.mean(predict(fit_boot, md_boot, type = "response"),
                      w = md_boot$weights)
  Odds1 <- P1 / (1 - P1)

  #Return marginal odds ratio
  return(Odds1 / Odds0)
}

boot_est <- boot(pair_ids, est_fun, R = 499)
boot_est
boot.ci(boot_est, type = "bca")

## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Cox Regression for marginal HR
coxph(Surv(Y_S) ~ A, data = md, robust = TRUE, 
      cluster = subclass)

## ---- eval = FALSE------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#  #Block bootstrap confidence interval
#  # library(boot)
#  
#  pair_ids <- levels(md$subclass)
#  
#  est_fun <- function(pairs, i) {
#  
#    #Compute number of times each pair is present
#    numreps <- table(pairs[i])
#  
#    #For each pair p, copy corresponding md row indices numreps[p] times
#    ids <- unlist(lapply(pair_ids[pair_ids %in% names(numreps)],
#                         function(p) rep(which(md$subclass == p),
#                                                numreps[p])))
#  
#    #Subset md with block bootstrapped ids
#    md_boot <- md[ids,]
#  
#    #Effect estimation
#    cox_fit_boot <- coxph(Surv(Y_S) ~ A, data = md_boot)
#  
#    #Compute the marginal HR by exponentiating the coefficient
#    #on treatment
#    HR <- exp(coef(cox_fit_boot)["A"])
#  
#    #Return the HR
#    return(HR)
#  }
#  
#  boot_est <- boot(pair_ids, est_fun, R = 499)
#  boot_est
#  boot.ci(boot_est, type = "bca")

## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
mNNr <- matchit(A ~ X1 + X2 + X3 + X4 + X5 + 
                 X6 + X7 + X8 + X9, data = d,
               link = "linear.logit", caliper = .1,
               ratio = 3, replace = TRUE)

mNNr

#match.data output
md <- match.data(mNNr)
nrow(md)
head(md)

#get_matches output
gm <- get_matches(mNNr)
nrow(gm)
head(gm)


## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Number of time control units are rematched
table(table(gm$id[gm$A == 0]))

## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Number of control units in each match stratum
table(table(gm$subclass[gm$A == 0]))

## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#match.data() output
fit1md <- lm(Y_C ~ A, data = md, weights = weights)

coeftest(fit1md, vcov. = vcovHC)

## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#get_matches() output
fit1gm <- lm(Y_C ~ A, data = gm, weights = weights)

coeftest(fit1gm, vcov. = vcovCL, cluster = ~subclass + id)

## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Full bootstrap confidence interval
# library(boot)

est_fun <- function(data, i) {
  #Matching function
  mNNr_boot <- matchit(A ~ X1 + X2 + X3 + X4 + X5 + 
                         X6 + X7 + X8 + X9, data = data[i,],
                       link = "linear.logit", caliper = .1,
                       ratio = 3, replace = TRUE)
  md_boot <- match.data(mNNr_boot)
  
  #Effect estimation
  fit_boot <- lm(Y_C ~ A, data = md_boot, weights = weights)
  
  #Return the coefficient on treatment
  return(coef(fit_boot)["A"])
}

boot_est <- boot(d, est_fun, R = 499)
boot_est
boot.ci(boot_est, type = "perc")

## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
fit2md <- lm(Y_C ~ A + X1 + X2 + X3 + X4 + X5 + 
                 X6 + X7 + X8 + X9, data = gm, 
             weights = weights)

coeftest(fit1gm, vcov. = vcovCL, cluster = ~subclass + id)["A",,drop = FALSE]

## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
fit3gm <- glm(Y_B ~ A, data = gm, weights = weights,
              family = quasibinomial(link = "logit"))

coeftest(fit3gm, vcov. = vcovCL, cluster = ~ subclass + id)
exp(coef(fit3gm)) #OR

## ---- eval = FALSE------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#  #Bootstrap confidence intervals
#  # library(boot)
#  
#  est_fun <- function(data, i) {
#    #Matching function
#    mNNr_boot <- matchit(A ~ X1 + X2 + X3 + X4 + X5 +
#                   X6 + X7 + X8 + X9, data = data[i,],
#                 link = "linear.logit", caliper = .1,
#                 ratio = 3, replace = TRUE)
#    md_boot <- match.data(mNNr_boot)
#  
#    #Fitting the model
#    fit_boot <- glm(Y_B ~ A + X1 + X2 + X3 + X4 + X5 +
#                      X6 + X7 + X8 + X9, data = md_boot,
#                    family = quasibinomial(link = "logit"),
#                    weights = weights)
#  
#    #Estimate potential outcomes for each unit
#    md_boot$A <- 0
#    P0 <- weighted.mean(predict(fit_boot, md_boot, type = "response"),
#                        w = md_boot$weights)
#    Odds0 <- P0 / (1 - P0)
#  
#    md_boot$A <- 1
#    P1 <- weighted.mean(predict(fit_boot, md_boot, type = "response"),
#                        w = md_boot$weights)
#    Odds1 <- P1 / (1 - P1)
#  
#    #Return marginal odds ratio
#    return(Odds1 / Odds0)
#  }
#  
#  boot_est <- boot(d, est_fun, R = 4999)
#  boot_est
#  boot.ci(boot_est, type = "bca")

## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Austin & Cafri's (2020) SE estimator
fs <- coxph(Surv(Y_S) ~ A, data = gm, robust = TRUE, 
      weights = weights, cluster = subclass)
Vs <- fs$var
ks <- nlevels(gm$subclass)

fi <- coxph(Surv(Y_S) ~ A, data = gm, robust = TRUE, 
      weights = weights, cluster = id)
Vi <- fi$var
ki <- length(unique(gm$id))

fc <- coxph(Surv(Y_S) ~ A, data = gm, robust = TRUE, 
      weights = weights)
Vc <- fc$var
kc <- nrow(gm)

#Compute the variance
V <- (ks/(ks-1))*Vs + (ki/(ki-1))*Vi - (kc/(kc-1))*Vc

#Sneak it back into the fit object
fc$var <- V

fc

## ---- include=FALSE-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#In case optmatch goes offline, don't run lines below
if (!requireNamespace("optmatch", quietly = TRUE)) knitr::opts_chunk$set(eval = FALSE)

## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
mF <- matchit(A ~ X1 + X2 + X3 + X4 + X5 + 
                 X6 + X7 + X8 + X9, data = d,
              method = "full", estimand = "ATE")

mF

md <- match.data(mF)

head(md)

## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
fit1 <- lm(Y_C ~ A, data = md, weights = weights)

coeftest(fit1, vcov. = vcovCL, cluster = ~subclass)

## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Estimating a covariate-adjusted marginal effect
#with treatment-covariate interactions

#Create a new dataset for centered variables
md_cen <- md

covs_to_center <- c("X1", "X2", "X3", "X4", "X5",
                    "X6", "X7", "X8", "X9")
md_cen[covs_to_center] <- scale(md_cen[covs_to_center], 
                                scale = FALSE)

#Fit the model with every covariate interacting with treatment
fit2 <- lm(Y_C ~ A * (X1 + X2 + X3 + X4 + X5 + 
                      X6 + X7 + X8 + X9),
           data = md_cen, weights = weights)

#Only output the intercept and coefficient on treatment
coeftest(fit2, vcov. = vcovCL, cluster = ~subclass)[1:2,]

## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
fit3 <- glm(Y_B ~ A, data = md, weights = weights,
            family = quasibinomial(link = "logit"))

coeftest(fit3, vcov. = vcovCL, cluster = ~subclass)
exp(coef(fit3)) #OR

## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
coxph(Surv(Y_S) ~ A, data = md, robust = TRUE, 
      weights = weights, cluster = subclass)

## ---- include=FALSE, eval=TRUE------------------------------------------------------------------------------------------------------------------------------------------------------------------------
knitr::opts_chunk$set(eval = TRUE)

## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
mS <- matchit(A ~ X1 + X2 + X3 + X4 + X5 + 
                 X6 + X7 + X8 + X9, data = d,
              method = "subclass", estimand = "ATT",
              subclass = 8)

mS

md <- match.data(mS)

head(md)

## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
fit1 <- lm(Y_C ~ A, data = md, weights = weights)

coeftest(fit1, vcov. = vcovHC)

## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
fit2 <- lm(Y_C ~ subclass + subclass:A - 1, data = md)

#Within-subclass effects
# coeftest(fit2, vcov. = vcovHC)

## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Subclass weights for ATT
sub_w <- with(md, c(rep(0, nlevels(subclass)), 
                    table(subclass[A==1])/sum(A==1)))

#Subclass weights for ATE (requires estimand = "ATE" in matchit())
# sub_w <- with(md, c(rep(0, nlevels(subclass)), 
#                     table(subclass)/nrow(md)))

#Marginal effect
(est <- weighted.mean(coef(fit2), sub_w))

#SE of marginal effect
(se <- sqrt(drop(sub_w %*% vcovHC(fit2) %*% sub_w)))

#CI
c(ci_low = est - 1.96*se, ci_hi = est + 1.96*se)

## ----eval=FALSE---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#  #Using margins() from margins
#  summary(margins::margins(fit2, variables = "A",
#                           data = md[md$A == 1,],
#                           vcov = vcovHC(fit2)))
#  #For ATE, omit the second line.

## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
fit3 <- glm(Y_B ~ A, data = md, weights = weights,
            family = quasibinomial(link = "logit"))

coeftest(fit3, vcov. = vcovHC)
exp(coef(fit3))

## ---- eval = FALSE------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#  #Bootstrap confidence intervals
#  library(boot)
#  
#  est_fun <- function(data, i) {
#    #Subclassification function
#    mS_boot <- matchit(A ~ X1 + X2 + X3 + X4 + X5 +
#                         X6 + X7 + X8 + X9, data = data[i,],
#                       method = "subclass", estimand = "ATT",
#                       subclass = 8)
#    md_boot <- match.data(mS_boot)
#  
#    #Fitting the model
#    fit_boot <- glm(Y_B ~ A * (subclass + X1 + X2 + X3 + X4 + X5 +
#                      X6 + X7 + X8 + X9), data = md_boot,
#                    family = quasibinomial(link = "logit"))
#  
#    #Estimate potential outcomes for each unit
#  
#    ## Subset to just the treated for the ATT; remove this for the ATE
#    md_boot <- md_boot[md_boot$A == 1,]
#    ##
#  
#    md_boot$A <- 0
#    P0 <- mean(predict(fit_boot, md_boot, type = "response"))
#    Odds0 <- P0 / (1 - P0)
#  
#    md_boot$A <- 1
#    P1 <- mean(predict(fit_boot, md_boot, type = "response"))
#    Odds1 <- P1 / (1 - P1)
#  
#    #Return marginal odds ratio
#    return(Odds1 / Odds0)
#  }
#  
#  boot_est <- boot(d, est_fun, R = 4999)
#  boot_est
#  boot.ci(boot_est, type = "bca")

## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
coxph(Surv(Y_S) ~ A, data = md, robust = TRUE, 
      weights = weights)

## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
mF <- matchit(A ~ X1 + X2 + X3 + X4 + X5 + 
                 X6 + X7 + X8 + X9, data = d,
              method = "full", estimand = "ATT")

md <- match.data(mF)

fitcond <- lm(Y_B ~ A + X1 + X2 + X3 + X4 + X5 + 
                 X6 + X7 + X8 + X9, data = md, 
             weights = weights)

coeftest(fitcond, vcov. = vcovCL, cluster = ~subclass)["A",,drop = FALSE]
exp(coef(fitcond)["A"])

## ---- eval = FALSE------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#  #Generatng data similar to Austin (2009) for demonstrating treatment effect estimation
#  gen_X <- function(n) {
#    X <- matrix(rnorm(9 * n), nrow = n, ncol = 9)
#    X[,5] <- as.numeric(X[,5] < .5)
#    X
#  }
#  
#  #~20% treated
#  gen_A <- function(X) {
#    LP_A <- - 1.2 + log(2)*X[,1] - log(1.5)*X[,2] + log(2)*X[,4] - log(2.4)*X[,5] + log(2)*X[,7] - log(1.5)*X[,8]
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
#  # Binary outcome
#  gen_Y_B <- function(A, X) {
#    LP_B <- -2 + log(2.4)*A + log(2)*X[,1] + log(2)*X[,2] + log(2)*X[,3] + log(1.5)*X[,4] + log(2.4)*X[,5] + log(1.5)*X[,6]
#    P_B <- plogis(LP_B)
#    rbinom(length(A), 1, P_B)
#  }
#  #Conditional:
#  #  OR:   2.4
#  #  logOR: .875
#  #Marginal:
#  #  RD:    .144
#  #  RR:   1.54
#  #  logRR: .433
#  #  OR:   1.92
#  #  logOR  .655
#  
#  # Survival outcome
#  gen_Y_S <- function(A, X) {
#    LP_S <- -2 + log(2.4)*A + log(2)*X[,1] + log(2)*X[,2] + log(2)*X[,3] + log(1.5)*X[,4] + log(2.4)*X[,5] + log(1.5)*X[,6]
#    sqrt(-log(runif(length(A)))*2e4*exp(-LP_S))
#  }
#  #Conditional:
#  #  HR:   2.4
#  #  logHR: .875
#  #Marginal:
#  #  HR:   1.57
#  #  logHR: .452
#  
#  set.seed(19599)
#  
#  n <- 2000
#  X <- gen_X(n)
#  A <- gen_A(X)
#  
#  Y_C <- gen_Y_C(A, X)
#  Y_B <- gen_Y_B(A, X)
#  Y_S <- gen_Y_S(A, X)
#  
#  d <- data.frame(A, X, Y_C, Y_B, Y_S)

