## Splitting functions for survival trees   

##  Evaluation
.surveval <- function(y, x, wt, params){
  tfit <- coxph(Surv(y[, 1], y[, 2]) ~ 1)
  wEvents <- tfit$nevent
  dev0 <- (-2)*tfit$loglik
  list(label = wEvents, deviance = dev0) }

##	Initialize fonction
.survinit <- function(y, params, wt){
  sfun <- function(yval, dev, wt, ylevel, digits){
    paste("events = ", round(yval),
          ", deviance = ", format(signif(dev, digits)),
          sep = '')}
  
  ttfun <- function (yval, dev, wt, ylevel, digits, n, use.n = TRUE) {
    # paste(yval,"\n", n, sep = '')
    paste(yval, n, sep = '/')
  }
  
  environment(sfun) <- .GlobalEnv
  environment(ttfun) <- .GlobalEnv
  list(y = cbind(y[, 1], y[, 2]), params = 0, numresp = 1, numy = 2, summary = sfun, text = ttfun)
}

##	Spliting function methode base (M0 - Logrank simple)
.survsplitLR <- function(y, wt, x, params, continuous){
  if(continuous)
  {
    #	continuous x variable: do all the logistic regressions
    n <- length(y[, 1])
    goodness <- double(n-1)
    direction <- goodness
    temp <- rep(0, n)
    for(i in 1:(n-1))
    {
      temp[i] <- 1
      if(x[i] != x[i+1])
      {
        tfit <- coxph(Surv(y[, 1], y[, 2]) ~ strata(log(wt)) + temp, robust = FALSE)
#         nbfit <- summary(tfit)$nevent
#         ngrp <- length(tfit$coef)
        goodness[i] <- tfit$score
        direction[i] <- sign(tfit$coef[1])
      }
    }
  }
  else
  {
    #	Categorical variable
    #	First, find out what order to put the categories in, which
    #	will be the order of the coefficients in this model
    x = x[, drop = TRUE]
    tfit <- glm(y[, 2] ~ factor(x), family = 'binomial')
    ngrp <- length(tfit$coef)
    direction <- rank(rank(tfit$coef) + runif(ngrp, 0, 0.1)) # break ties
    xx <- direction[match(x, sort(unique(x)))] # relabel from small to large
    goodness <- double(length(direction) - 1)
    for(i in 1:length(goodness))
    {
      tfit <- coxph(Surv(y[, 1], y[, 2]) ~ strata(log(wt)) + I(xx > i), robust = FALSE)
#       nbfit <- summary(tfit)$nevent
      goodness[i] <- tfit$score
      
    }
  }
  
  list(goodness = goodness, direction = direction)
}

##::##	Spliting function methode cure-rate (M1)
.survsplitR2 <- function(y, wt, x, params, continuous){
  if(continuous)
  {
    #	continuous x variable: do all the Cox regressions
    n <- length(y[, 1])
    goodness <- double(n-1)
    direction <- goodness
    temp <- rep(0, n)
    for(i in 1:(n-1))
    {
      temp[i] <- 1
      if(x[i] != x[i+1])
      {
        tfit <- coxph(Surv(y[, 1], y[, 2]) ~ temp, robust = T)

        goodness[i] <- PseudoR2.Cure(ygene = temp, ydelai = y[, 1], yetat = y[, 2], strate = log(wt))
        direction[i] <- sign(tfit$coef[1])
      }
    }
  }else
  {
    #	Categorical variable
    #	First, find out what order to put the categories in, which

    x = x[, drop = TRUE]
    tfit <- glm(y[, 2] ~ factor(x), family = 'binomial')
    ngrp <- length(tfit$coef)
    direction <- rank(rank(tfit$coef) + runif(ngrp, 0, 0.1)) # break ties
    xx <- direction[match(x, sort(unique(x)))] # relabel from small to large
    goodness <- double(length(direction) - 1)
    for(i in 1:length(goodness))
    {

      goodness[i] <- PseudoR2.Cure(ygene = I(xx > i), ydelai = y[, 1], yetat = y[, 2], strate = log(wt))
      
    }
  }
  
  list(goodness = goodness, direction = direction)
}
