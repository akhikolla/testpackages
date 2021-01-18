# Start bgtCI() functions
###############################################################################

##################################################################
# bgtAC() function                                               #
##################################################################

"bgtAC" <-
  function(n, y, s, conf.level = 0.95, alternative = "two.sided") {
    alpha <- 1 - conf.level
    
    est.int <- (y + (qnorm(1 - alpha/2)^2)/2)/(n + (qnorm(1 - alpha/2))^2)
    est.int1s <- (y + (qnorm(1 - alpha)^2)/2)/(n + (qnorm(1 - alpha))^2)
    
    
    if (alternative == "two.sided") {
      AC.se <- (qnorm(1 - alpha/2))*sqrt((est.int*(1 - est.int))/(n + (qnorm(1 - alpha/2))^2))
      KI.int.l <- est.int - AC.se
      KI.int.u <- est.int + AC.se
      if (KI.int.u > 1) {KI.int.u <- 1}
      if (KI.int.l < 0) {KI.int.l <- 0}
      CI <- c(1 - (1 - KI.int.l)^(1/s), 1 - (1 - KI.int.u)^(1/s))
    } 
    
    else{if (alternative == "less") {
      AC.se <- (qnorm(1 - alpha))*sqrt((est.int1s*(1 - est.int1s))/(n + (qnorm(1 - alpha))^2))
      KI.int.u <- est.int1s + AC.se
      if (KI.int.u > 1) {KI.int.u <- 1}
      CI <- c(0, 1 - (1 - KI.int.u)^(1/s))
    } 
      
      else{if (alternative == "greater") {
        AC.se <- (qnorm(1 - alpha))*sqrt((est.int1s*(1 - est.int1s))/(n + (qnorm(1 - alpha))^2))
        KI.int.l <- est.int1s - AC.se
        if (KI.int.l < 0) {KI.int.l <- 0}
        CI <- c(1 - (1 - KI.int.l)^(1/s),1)
      }
        
        else{stop("argument alternative misspecified")}}}
    
    CI
  }

##################################################################
# bgtBlaker() function                                           #
##################################################################

"bgtBlaker" <-
  function(n, y, s, conf.level = 0.95, alternative = "two.sided")
  {
    
    # # # from the S code given in Blaker(2000), slightly changed # # #
    
    tolerance <- 1e-04
    
    acceptbin <- function(y,n,p)
    {
      p1 <- 1 - pbinom(y - 1, n, p)
      p2 <- pbinom(y, n, p)
      a1 <- p1 + pbinom( qbinom(p1,n,p) - 1, n, p )
      a2 <- p2 + 1 - pbinom( qbinom(1 - p2, n, p), n, p )
      return(min(a1,a2))
    }
    
    lower <- 0
    upper <- 1
    
    if (y != 0)
    {lower <- qbeta((1 - conf.level)/2, y, n - y + 1)
    {while (acceptbin(y, n, lower + tolerance) < (1 - conf.level))
      lower <- lower + tolerance}
    }
    
    if (y != n)
    {upper <- qbeta(1 - (1 - conf.level)/2, y + 1, n - y)
    {while (acceptbin(y, n, upper - tolerance) < (1 - conf.level))
      upper <- upper - tolerance}
    }
    CI <- c(1 - (1 - lower)^(1/s), 1 - (1 - upper)^(1/s))
    
    CI  
  }




##################################################################
# bgtCP() function                                               #
##################################################################

"bgtCP" <-
  function(n, y, s, conf.level = 0.95, alternative = "two.sided")
  {
    lower <- 0
    upper <- 1
    if (alternative == "two.sided")
    {
      if (y != 0)
      {lower <- qbeta((1 - conf.level)/2, y, n - y + 1)}
      
      if (y != n)
      {upper <- qbeta(1 - (1 - conf.level)/2, y + 1, n - y)}
    }
    
    if (alternative == "less")
    {
      if (y != n)
      {upper <- qbeta(1 - (1 - conf.level), y + 1, n - y)}
    }
    
    if (alternative == "greater")
    {
      if (y != 0)
      {lower <- qbeta((1 - conf.level), y, n - y + 1)}
    }
    
    estimate <- 1 - (1 - y/n)^(1/s)
    
    CI <- c(1 - (1 - lower)^(1/s),1 - (1 - upper)^(1/s))
    
    CI   
  }




##################################################################
# bgtSOC() function                                              #
##################################################################

"bgtSOC" <-
  function(n, s, y, conf.level = 0.95, alternative = "two.sided")
    
  {
    esti <- y/n
    kappa <- qnorm(conf.level)
    eta <- (kappa^2)/3 + 1/6
    gamma1 <- ((13/18)*kappa^2 + 17/18)*(-1)
    gamma2 <- (kappa^2)/18 + 7/36
    
    midpo <- (y + eta)/(n + 2*eta)
    
    if (alternative == "less")
    {upper <- midpo + kappa * sqrt(esti*(1 - esti) + (gamma1*esti*(1 - esti) + gamma2)/n)/sqrt(n)
    CI <- c(0 ,upper)
    if (y == n || upper > 1) {CI <- c(0,1)}
    else{CI <- c(0 ,upper)}
    }
    
    if (alternative == "greater")
    {CI <- c(midpo - kappa*sqrt(esti*(1 - esti) + (gamma1*esti*(1 - esti) + gamma2)/n)/sqrt(n), 1)
    if (y == 0) {CI <- c(0,1)} }
    
    if (alternative == "two.sided")
    {
      kappa <- qnorm(1 - (1 - conf.level)/2)
      eta <- (kappa^2)/3 + 1/6
      gamma1 <- ((13/18)*kappa^2 + 17/18)*(-1)
      gamma2 <- (kappa^2)/18 + 7/36
      
      lower <- midpo - kappa*sqrt(esti*(1 - esti) + (gamma1*esti*(1 - esti) + gamma2)/n)/sqrt(n)  
      upper <- midpo + kappa*sqrt(esti*(1 - esti) + (gamma1*esti*(1 - esti) + gamma2)/n)/sqrt(n)
      
      if (y == 0) {CI <- c(0,upper)} 
      else{if (y == n || upper > 1) {CI <- c(lower,1)}
        else{CI <- c(lower, upper)}}
    }
    
    CI2 <- c(1 - (1 - CI[1])^(1/s),  1 - (1 - CI[2])^(1/s))
    CI2
  }




##################################################################
# bgtWald() function                                             #
##################################################################

"bgtWald" <-
  function(n, y, s, conf.level = 0.95, alternative = "two.sided")
    
  {
    if (y > n) {stop("number of positive tests y can not be greater than number of tests n")}
    th <- y/n
    esti <- 1 - (1 - th)^(1/s)
    var.esti <- (1 - (1 - esti)^s)/(n*(s^2)*(1 - esti)^(s - 2))
    alpha <- 1 - conf.level
    
    if (alternative == "two.sided") {
      snquant <- qnorm(p = 1 - alpha/2, mean = 0, sd = 1, lower.tail = TRUE)
      CI <- c(esti - snquant*sqrt(var.esti), esti + snquant*sqrt(var.esti))
    }
    else{if (alternative == "less") {
      snquant <- qnorm(p = 1 - alpha, mean = 0, sd = 1, lower.tail = TRUE)
      CI <- c(0, esti + snquant*sqrt(var.esti))
    }
      else {if (alternative == "greater") {
        snquant <- qnorm(p = 1 - alpha, mean = 0, sd = 1, lower.tail = TRUE)
        CI <- c(esti - snquant*sqrt(var.esti), 1)
      }
        else {stop("argument alternative mis-specified")}}}
    CI
  }

##################################################################
# bgtWilson() function                                           #
##################################################################

"bgtWilson" <-
  function(n, y, s, conf.level = 0.95, alternative = "two.sided")
  { 
    alpha <- 1 - conf.level 
    th <- y/n
    est.int <- (y + (qnorm(1 - alpha/2)^2)/2)/(n + (qnorm(1 - alpha/2))^2)
    est.int1s <- (y + (qnorm(1 - alpha)^2)/2)/(n + (qnorm(1 - alpha))^2)
    
    if (alternative == "two.sided") {
      w.se <- ((qnorm(1 - alpha/2))*sqrt(n*th*(1 - th) + (qnorm(1 - alpha/2)^2)/4))/(n + qnorm(1 - alpha/2)^2)
      KI.int.l <- est.int - w.se
      KI.int.u <- est.int + w.se
      if (KI.int.u > 1) {KI.int.u <- 1}
      if (KI.int.l < 0) {KI.int.l <- 0}
      KI <- c(1 - (1 - KI.int.l)^(1/s), 1 - (1 - KI.int.u)^(1/s))
    }
    
    else{if (alternative == "less") {
      w.se <- ((qnorm(1 - alpha))*sqrt(n*th*(1 - th) + (qnorm(1 - alpha)^2)/4))/(n + qnorm(1 - alpha)^2)
      KI.int.u <- est.int1s + w.se
      if (KI.int.u > 1) {KI.int.u <- 1}
      KI <- c( 0, 1 - (1 - KI.int.u)^(1/s) )
    }
      
      else{if (alternative == "greater") {
        w.se <- ((qnorm(1 - alpha))*sqrt(n*th*(1 - th) + (qnorm(1 - alpha)^2)/4))/(n + qnorm(1 - alpha)^2)
        KI.int.l <- est.int1s - w.se
        if (KI.int.l < 0) {KI.int.l <- 0}
        KI <- c(1 - (1 - KI.int.l)^(1/s), 1)
      }
        
        else{stop("argument alternative misspecified")}}}
    
    KI
  }




##################################################################
# bgtCI() function                                               #
##################################################################

# Brianna Hitt - 02/10/2020
# Minor changes were made to the capitalization of confidence interval 
#   method names - only the first letter is capitalized for methods 
#   named after people; otherwise, all lowercase

# Brianna Hitt - 04.02.2020
# Changed cat()/print() to warning()

"bgtCI" <-
  function(n, s, y, conf.level = 0.95, alternative = "two.sided", method = "CP")
    
  {
    if(length(n)!=1 || (n<1 | abs(round(n)-n) > 1e-07)){stop("number of groups n must be specified as a single integer > 0")}
    if(length(s)!=1 || (s<1 | abs(round(s)-s) > 1e-07)){stop("group size s must be specified as a single integer > 0")}
    if(length(s)!=1 || (y<0 | abs(round(y)-y) > 1e-07)){stop("observed number of positive groups y must be specified as a single integer>0")}
    if(y>n) {stop("number of positive tests y can not be greater than number of groups n")}
    if(length(conf.level)!=1 || conf.level<0 || conf.level>1){stop("conf.level must be a positive number between 0 and 1")}
    
    method<-match.arg(method, choices=c("CP","Blaker","AC","score","Wald","soc"))
    alternative<-match.arg(alternative, choices=c("two.sided","less","greater"))
    
    estimate=1-(1-y/n)^(1/s)
    
    switch(method,
           "CP"={conf.int<-bgtCP(n=n, s=s, y=y, conf.level=conf.level, alternative=alternative)},
           
           "Blaker"={
             if(alternative=="less" || alternative=="greater"){
               warning("The Blaker CI is inherently two.sided")
               conf.int<-c(NA, NA)}
             if(alternative=="two.sided")
             {conf.int<-bgtBlaker(n=n, s=s, y=y, conf.level=conf.level)}
           },
           
           "score"={conf.int<-bgtWilson(n=n, s=s, y=y, conf.level=conf.level, alternative=alternative)},
           
           "AC"={conf.int<-bgtAC(n=n, s=s, y=y, conf.level=conf.level, alternative=alternative)},
           
           "Wald"={conf.int<-bgtWald(n=n, s=s, y=y, conf.level=conf.level, alternative=alternative)},
           
           "soc"={conf.int<-bgtSOC(n=n, s=s, y=y, conf.level=conf.level, alternative=alternative)}
    )
    
    out<-list(conf.int=conf.int,
              estimate=estimate,
              method=method,
              conf.level=conf.level,
              alternative=alternative)
    class(out)<-"bgtCI"
    out
  }




# Start bgtvs() functions
###############################################################################

##################################################################
# bgtvs() function                                               #
##################################################################

"bgtvs" <- function
(n, s, y, conf.level = 0.95, alternative = "two.sided", maxiter = 100)
{
  
  alternative <- match.arg(alternative, choices = c("two.sided", "less", "greater"))
  
  # error checking:
  
  if (length(y) != length(n) || length(n) != length(s))
  {stop("vectors n, s, y must have exactly the same length")}
  
  if ( any(y > n))
  {stop("values of y must be smaller than or equal to the corresponding values n") }
  
  if (any(abs(round(n) - n) > 1e-07) | any(n < 1))
  {stop("number of groups n must be a vector of integer values > 0")}
  
  if (any(abs(round(s) - s) > 1e-07) | any(s < 1))
  {stop("group sizes s must be a vector of integer values > 0")}
  
  if (any(abs(round(y) - y) > 1e-07) | any(y < 0))
  {stop("number of positive groups y must be a vector of non-negative integer values > 0")}
  
  
  if (conf.level <= 0 | conf.level >= 1)
  {stop("conf.level must be a numeric value between 0, and 1, usually e.g. 0.95")}
  
  
  m <- length(y)
  alpha <- 1 - conf.level
  
  # # # Likelihood according to Hepworth 1996 (1):
  
  likHepI <- function(n, s, y, p)
  {
    prod(choose(n,y) * ((1 - (1 - p)^s )^y) * ((1 - p)^(s*(n - y))))
  }
  
  # function to estimate the MLE:
  
  estBGTvgs2 <- function(n, s, y)
  {
    if (length(y) != length(n) || length(n) != length(s)) 
    {stop("vectors n, s, y must have exactly the same length")}
    
    # total number of individuals in the trial
    
    # the current iteration method is sensitive 
    # to give wrong results for too small tolerances
    
    tol <- 1e-10	
    total <- sum(n*s)
    
    # Hepworth(1996), equation 2
    
    mindiff <- function(n, s, y, p, total)
    {total - sum((s*y) / (1 - (1 - p)^(s)))}
    
    # check the case all y=0 (equa 2 not defined),
    # else iterate p
    
    if (all(y == 0))
    {esti <- 0}
    
    else
    {
      maxiter <- 100
      dir <- 1
      crit <- numeric(length = maxiter)
      
      if (any(y == 0))
      {esti <- 0 + tol}
      else
      {esti <- 0}
      
      crit[1] <- mindiff(n = n, s = s, y = y, p = esti, total = total)
      step <- 0.5
      
      for (i in 2:maxiter)
      {
        crit[i] <- mindiff(n = n, s = s, y = y, p = esti, total = total)
        if (sign(crit[i - 1]) != sign(crit[i]))
        {dir <- dir*(-1); step <- step/2}
        if (esti > 1 | esti < 0)
        {dir <- dir*(-1)}
        esti <- esti + dir*step
      }
    }
    return(esti)
  }
  
  # # # Point estimate for the observed outcome:
  # Pobs (Hepworth)
  
  point.esti <- estBGTvgs2(n = n, s = s, y = y)
  
  
  # # # calculate the MLE for all possible combinations of y1, ..., ym
  
  possibley <- list()
  
  for (i in 1:length(n))
  {
    possibley[[i]] <- 0:n[i]
  }
  allComb <- expand.grid(possibley)
  
  
  MLEComb <- numeric(length = nrow(allComb))
  
  for (comb in 1:nrow(allComb))
  {ytemp <- as.numeric(allComb[comb,])
  MLEComb[comb] <- estBGTvgs2(n = n, s = s, y = ytemp)
  }
  
  out <- cbind(allComb,MLEComb)
  
  # # # order the event space by their associated MLE
  
  outlower <- as.matrix( out[which(out$MLEComb >= point.esti), 1:m] )
  
  outupper <- as.matrix( out[which(out$MLEComb <= point.esti), 1:m])
  
  # Iteration functions for the upper and lower bound:
  
  itupperlimit <- function(pstart, outupper, alpha, n, s, maxiter = maxiter)
  {
    dir <- 1
    step <- 0.1
    
    crit <- numeric(length = maxiter)
    
    crit[1] <- alpha - sum(apply(X = outupper, MARGIN = 1, FUN = likHepI, 
                                 n = n, s = s, p = pstart))
    
    p <- pstart + dir*step
    
    for (i in 2:maxiter)
    {
      crit[i] <- alpha - sum(apply(X = outupper, MARGIN = 1, FUN = likHepI, 
                                   n = n, s = s, p = p))
      if (sign(crit[i - 1]) != sign(crit[i]))
      {dir <- dir*(-1); step <- step/2}
      
      p <- min(p + dir*step, 1)
    }
    
    return(pu = p)
  }
  
  
  itlowerlimit <- function(pstart, outlower, alpha, n, s, maxiter = maxiter)
  {
    dir <- (-1)
    step <- 0.1
    
    crit <- numeric(length = maxiter)
    
    crit[1] <- alpha - sum(apply(X = outlower, MARGIN = 1, FUN = likHepI, 
                                 n = n, s = s, p = pstart))
    
    p <- pstart + dir*step
    
    for (i in 2:maxiter)
    {
      crit[i] <- alpha - sum(apply(X = outlower, MARGIN = 1, FUN = likHepI, n = n, s = s, p = p))
      if (sign(crit[i - 1]) != sign(crit[i]))
      {dir <- dir*(-1); step <- step/2}
      
      p <- max(p + dir*step,0)
    }
    
    return(pl = p)
  }
  
  # # # Calculation of Bounds
  
  if (alternative == "two.sided")
  {
    upper <- itupperlimit(pstart = point.esti, outupper = outupper, alpha = alpha/2, n = n, s = s, maxiter = 100)
    lower <- itlowerlimit(pstart = point.esti, outlower = outlower, alpha = alpha/2, n = n, s = s, maxiter = 100)
  }
  
  if (alternative == "less")
  {
    upper <- itupperlimit(pstart = point.esti, outupper = outupper, alpha = alpha, n = n, s = s, maxiter = 100)
    lower <- 0
  }
  
  if (alternative == "greater")
  {
    upper = 1
    lower <- itlowerlimit(pstart = point.esti, outlower = outlower, alpha = alpha/2, n = n, s = s, maxiter = 100)
  }
  
  out <- list(conf.int = c(lower, upper),
              estimate = point.esti, 
              conf.level = conf.level,
              alternative = alternative,
              input = rbind("number of groups" = n,
                            "group size" = s,
                            "number of positive groups" = y)
  )
  
  class(out) <- "bgtvs"
  
  return(out)
}




##################################################################
# pooledBin() function                                           #
##################################################################

# Brad Biggerstaff
# Division of Vector-Borne Infectious Diseases
# National Center for Infectious Diseases
# Centers for Disease Control and Prevention
# P.O. Box 2087, Fort Collins, CO  80522-2087
# (970) 221-6473 ... BBiggerstaff@cdc.gov
#
# In all of these functions:
#    m:  a vector of pool sizes
#    n:  a vector of the corresponding number of pools of sizes m
#    x:  a vector of the corresponding number of the n pools of size m that are positive
#    p:  the proportion
#
# Options include:
#   pt.est: specify the point estimate to compute, including the MLE ("mle"),
#           bias-corrected MLE ("bc-mle") [default], and MIR ("mir").
#       ci: specify the confidence interval to compute, including the score ("score"), s
#           kewness-corrected ("skew-score") [default], score ("score"), bias- and skewness-corrected
#           score ("bc-skew-score"), likelihood ratio test ("lrt"), Wald ("Wald"), and MIR ("mir")
#    alpha: 1-confidence level [default = 0.05]
#      tol: the tolerance for convergence to use in the algorithms
#
# Various confidence intervals are given, as noted by *.ci in the function name
#
# References:
#	Walter SD, Hildreth SW, Beaty BJ: Estimation of infection rates in population of
#   organisms using pools of variable size. Am J Epidemiol 1980, 112(1):124-128.
#
# Hepworth G: Estimation of proportions by group testing. PhD Dissertation.
#   Melbourne, Australia: The University of Melbourne; 1999.
#
# Biggerstaff BJ:  Confidence interval for the difference of proportions estmimated
#   from pooled samples.  JABES 2008, 13(4):478-496.
#
# Hepworth G, Biggerstaff BJ:  Bias correction in estimating proportions
#   by pooled testing.  JABES 2017, to appear.
#
#
################################################################################
# Toy examples
# ------------
#
#> x1 <- c(1,0,0,0,0)
#> m1 <- c(10,4,1,25,50)
#> n1 <- c(5,1,1,30,20)
#> x2 <- c(2,0,1,0,0)
#> m2 <- c(10,4,1,25,50)
#> n2 <- c(5,1,1,30,20)
#--------------
#> pooledBin(x1,m1,n1)
#PointEst    Lower    Upper
#  0.0005   0.0000   0.0027
#
#--------------
#> pooledBin(x1,m1,n1,scale=1000)
# PointEst  Lower  Upper Scale
#   0.5497 0.0319 2.6524  1000
#
#--------------
#> summary(pooledBin(x1,m1,n1),scale=1000)
#Estimation of Binomial Proportion for Pooled Data
#
# PointEst  Lower  Upper Scale
#   0.5497 0.0319 2.6524  1000
#
#Point estimator: Firth's Correction
#CI method: Skew-Corrected score (Gart)
#
#Number of individuals: 1805
#Number of pools: 57
#Number of positive pools: 1
#
#--------------
#> pooledBin(x2,m2,n2,pt.method="mle",ci.method="lrt")
#PointEst    Lower    Upper
#  0.0017   0.0004   0.0043
#--------------
#> pooledBinDiff(x1,m1,x2,m2,n1,n2)
#PointEst    Lower    Upper
# -0.0011  -0.0040   0.0014
#
#--------------
#> summary(pooledBinDiff(x1,m1,x2,m2,n1,n1), scale=1000)
#Estimation of Difference of Binomial Proportions for Pooled Data
#
# PointEst   Lower  Upper Scale
#  -1.1036 -3.9916 1.3646  1000
#
#Point estimator: Firth's Correction
#CI method: Skew-Corrected score (Gart)
#
#             PointEst  Lower  Upper Scale Individuals Pools Positive Pools
#Population 1   0.5497 0.0319 2.6524  1000        1805    57              1
#Population 2   1.6533 0.4425 4.4240  1000        1805    57              3
#
#--------------
#> pooledBinDiff(x1,m1,x2,m2,n1,n2,pt.method="mle", ci.method="lrt")
#PointEst    Lower    Upper
# -0.0011  -0.0039   0.0012

#################################################################################
#################################################################################
#
# General functions to use to access the various methods
#
#################################################################################
#################################################################################

# Brianna Hitt - 02/10/2020
# Minor changes were made to the capitalization of confidence interval 
#   method names - only the first letter is capitalized for methods 
#   named after people; otherwise, all lowercase

"pooledBin" <-
  function(x,m,n=rep(1,length(x)),
           pt.method = c("Firth","Gart","bc-mle","mle","mir"),
           ci.method = c("skew-score","bc-skew-score","score","lrt","Wald","mir"),
           scale=1, alpha=0.05, tol=.Machine$double.eps^0.5)
  {
    call <- match.call()
    pt.method <- match.arg(pt.method)
    ci.method <- match.arg(ci.method)
    if (ci.method == "mir" | pt.method == "mir") {
      ci.method <- "mir"
      pt.method <- "mir"
    }
    if (pt.method == "Gart") pt.method <- "bc-mle" # backward compatability
    switch(pt.method,
           "Firth" = {p <- pooledbinom.firth(x,m,n,tol)},
           "mle" = {p <- pooledbinom.mle(x,m,n,tol)},
           "bc-mle" = {p <- pooledbinom.cmle(x,m,n,tol)},
           "mir" = {p <- pooledbinom.mir(x,m,n)}
    )
    if (p < 0 & pt.method == "bc-mle") {
      pt.method <- "mle"
      warning("Bias-correction results in negative point estimate; using MLE\n")
      p <- pooledbinom.mle(x,m,n,tol)
    }
    switch(ci.method,
           "skew-score" = {ci.p <- pooledbinom.cscore.ci(x,m,n,tol,alpha)[2:3]},
           "bc-skew-score" = {ci.p <- pooledbinom.bcscore.ci(x,m,n,tol,alpha)[2:3]},
           "score" = {ci.p <- pooledbinom.score.ci(x,m,n,tol,alpha)[2:3]},
           "lrt" = {ci.p <- pooledbinom.lrt.ci(x,m,n,tol,alpha)[2:3]},
           "Wald" = {ci.p <- pooledbinom.wald.ci(x,m,n,tol,alpha)[2:3]},
           "mir" = {ci.p <- pooledbinom.mir.ci(x,m,n,tol,alpha)[2:3]}
    )
    structure(list(p = p, lcl = ci.p[1], ucl = ci.p[2], 
                   pt.method = pt.method, ci.method = ci.method, alpha = alpha, 
                   x = x, m = m, n = n, scale = scale, call = call), 
              class = "poolbin")
  }




##################################################################
# propCI() function                                              #
##################################################################

#' @title Confidence intervals for one proportion in group testing
#' 
#' @description Calculates point estimates and confidence intervals for a 
#' single proportion with group testing data. Methods are available for groups 
#' of equal or different sizes. 
#' 
#' @param x integer specifying the number of positive groups when groups are 
#' of equal size, or a vector specifying the number of positive groups among 
#' the \kbd{n} groups tested when group sizes differ. If the latter, this 
#' vector must be of the same length as the \kbd{m} and \kbd{n} arguments.
#' @param m integer specifying the common size of groups when groups are of 
#' equal size, or a vector specifying the group sizes when group sizes differ. 
#' If the latter, this vector must be of the same length as the \kbd{x} and 
#' \kbd{n} arguments. 
#' @param n integer specifying the number of groups when these groups are of 
#' equal size, or a vector specifying the corresponding number of groups of 
#' the sizes \kbd{m} when group sizes differ. If the latter, this vector must 
#' be of the same length as the \kbd{x} and \kbd{m} arguments.
#' @param pt.method character string specifying the point 
#' estimate to compute. Options include \kbd{"Firth"} for the 
#' bias-preventative, \kbd{"Gart"} and \kbd{"bc-mle"} for the 
#' bias-corrected MLE (where the latter allows for backward compatibility), 
#' and \kbd{"mle"} for the MLE.
#' @param ci.method character string specifying the confidence 
#' interval to compute. Options include "AC" for the Agresti-Coull interval, 
#' "bc-skew-score" for the bias- and skewness-corrected interval, "Blaker" 
#' for the Blaker interval, "CP" for the Clopper-Pearson interval, "exact" 
#' for the exact interval as given by Hepworth (1996), "lrt" for the 
#' likelihood ratio test interval, "score" for the Wilson score interval, 
#' "skew-score" for the skewness-corrected interval, "soc" for the 
#' second-order corrected interval, and "Wald" for the Wald interval. Note 
#' that the Agresti-Coull, Blaker, Clopper-Pearson, and second-order corrected 
#' intervals can only be calculated when \kbd{x}, \kbd{m}, and \kbd{n} are 
#' given as integers (equal group size case). 
#' @param conf.level confidence level of the interval.
#' @param alternative character string defining the alternative 
#' hypothesis, either \kbd{"two.sided"}, \kbd{"less"}, or \kbd{"greater"}.
#' @param maxiter the maximum number of steps in the iteration of 
#' confidence limits, for use only with the \kbd{"exact"} method when 
#' group sizes differ.
#' @param tol the accuracy required for iterations in internal functions, 
#' for use with asymptotic intervals when group sizes differ only.
#' 
#' @details Confidence interval methods include the Agresti-Coull 
#' (\kbd{ci.method="AC"}), bias- and skewness-corrected 
#' (\kbd{ci.method="bc-skew-score"}), Blaker (\kbd{ci.method="Blaker"}), 
#' Clopper-Pearson (\kbd{ci.method="CP"}), exact (\kbd{ci.method="exact"}), 
#' likelihood ratio test (\kbd{ci.method="lrt"}), Wilson score 
#' (\kbd{ci.method="score"}), skewness-corrected 
#' (\kbd{ci.method="skew-score"}), second-order corrected 
#' (\kbd{ci.method="soc"}), and Wald 
#' (\kbd{ci.method="Wald"}) intervals. The Agresti-Coull, Blaker, 
#' Clopper-Pearson, and second-order corrected intervals are available 
#' only for the equal group size case.
#' 
#' Point estimates available include the MLE (\kbd{pt.method="mle"}), 
#' bias-corrected MLE (\kbd{pt.method="Gart"} or \kbd{pt.method="bc-mle"}), 
#' and bias-preventative (\kbd{pt.method="Firth"}). Only the MLE method 
#' is available when calculating the Clopper-Pearson, Blaker, Agresti-Coull, 
#' second-order corrected, or exact intervals.
#' 
#' \subsection{Equal group sizes}{
#' Computation of confidence intervals for group testing with equal group 
#' sizes are described in Tebbs & Bilder (2004) and Schaarschmidt (2007).} 
#' 
#' \subsection{Unequal group sizes}{
#' While the exact method is available when group sizes differ, 
#' the algorithm becomes computationally very expensive if the number of 
#' different groups, \kbd{n}, becomes larger than three. See Hepworth (1996) 
#' for additional details on the exact method and other methods for 
#' constructing confidence intervals in group testing situations. For 
#' computational details and simulation results of the remaining methods, 
#' see Biggerstaff (2008). See Hepworth & Biggerstaff (2017) for 
#' recommendations on the best point estimator methods.}
#'  
#' @return A list containing:
#' \item{conf.int}{a confidence interval for the proportion.}
#' \item{estimate}{the point estimator of the proportion.}
#' \item{pt.method}{the method used for point estimation.}
#' \item{ci.method}{the method used for confidence interval estimation.}
#' \item{conf.level}{the confidence level of the interval.}
#' \item{alternative}{the alternative specified by the user.}
#' \item{x}{the number of positive groups.}
#' \item{m}{the group sizes.}
#' \item{n}{the numbers of groups with corresponding group sizes \kbd{m}.}
#' 
#' @author This function is a combination of \code{bgtCI} and \code{bgtvs} 
#' written by Frank Schaarschmidt and \code{pooledBin} written by Brad 
#' Biggerstaff for the \code{binGroup} package. Minor modifications have been 
#' made for inclusion of the functions in the \code{binGroup2} package.
#' 
#' @references
#' \insertRef{Biggerstaff2008}{binGroup2}
#'   
#' \insertRef{Hepworth1996}{binGroup2}
#' 
#' \insertRef{Hepworth2017}{binGroup2}
#' 
#' \insertRef{Schaarschmidt2007}{binGroup2}
#' 
#' \insertRef{Tebbs2004}{binGroup2}
#' 
#' @seealso \code{\link{propDiffCI}} for confidence intervals for the 
#' difference of proportions in group testing, \code{\link{gtTest}} for 
#' hypothesis tests in group testing, \code{\link{gtPower}} for power 
#' calculations in group testing, and \code{\link{binom.test}} for an exact 
#' confidence interval and test.
#' 
#' @family estimation functions
#' 
#' @examples 
#' # Example from Tebbs and Bilder (2004):
#' #   3 groups out of 24 test positively; 
#' #   each group has a size of 7.
#' # Clopper-Pearson interval:
#' propCI(x=3, m=7, n=24, ci.method="CP", 
#'        conf.level=0.95, alternative="two.sided")
#'       
#' # Blaker interval:
#' propCI(x=3, m=7, n=24, ci.method="Blaker", 
#'        conf.level=0.95, alternative="two.sided")
#'       
#' # Wilson score interval: 
#' propCI(x=3, m=7, n=24, ci.method="score", 
#'        conf.level=0.95, alternative="two.sided")
#'       
#' # One-sided Clopper-Pearson interval:
#' propCI(x=3, m=7, n=24, ci.method="CP", 
#'        conf.level=0.95, alternative="less")
#' 
#' # Calculate confidence intervals with a group size of 1. 
#' #   These match those found using the binom.confint() 
#' #   function from the binom package.
#' propCI(x = 4, m = 1, n = 10, pt.method = "mle", 
#'        ci.method = "AC")
#' propCI(x = 4, m = 1, n = 10, pt.method = "mle", 
#'        ci.method = "score")
#' propCI(x = 4, m = 1, n = 10, pt.method = "mle", 
#'        ci.method = "Wald")
#' 
#' # Example from Hepworth (1996, table 5):
#' #   1 group out of 2 tests positively with 
#' #   groups of size 5; also, 
#' #   2 groups out of 3 test positively with 
#' #   groups of size 2.
#' propCI(x=c(1,2), m=c(5,2), n=c(2,3), ci.method="exact") 
#' 
#' # Recalculate the example given in
#' #   Hepworth (1996), table 5:
#' propCI(x=c(0,0), m=c(5,2), n=c(2,3), ci.method="exact")
#' propCI(x=c(0,1), m=c(5,2), n=c(2,3), ci.method="exact")
#' propCI(x=c(0,2), m=c(5,2), n=c(2,3), ci.method="exact")
#' propCI(x=c(0,3), m=c(5,2), n=c(2,3), ci.method="exact")
#' propCI(x=c(1,0), m=c(5,2), n=c(2,3), ci.method="exact")
#' propCI(x=c(1,1), m=c(5,2), n=c(2,3), ci.method="exact")
#' propCI(x=c(1,2), m=c(5,2), n=c(2,3), ci.method="exact")
#' propCI(x=c(1,3), m=c(5,2), n=c(2,3), ci.method="exact")
#' propCI(x=c(2,0), m=c(5,2), n=c(2,3), ci.method="exact")
#' propCI(x=c(2,1), m=c(5,2), n=c(2,3), ci.method="exact")
#' propCI(x=c(2,2), m=c(5,2), n=c(2,3), ci.method="exact")
#' propCI(x=c(2,3), m=c(5,2), n=c(2,3), ci.method="exact")
#' 
#' # Example with multiple groups of various sizes: 
#' #   0 out of 5 groups test positively with 
#' #   groups of size 1 (individual testing);
#' #   0 out of 5 groups test positively with 
#' #   groups of size 5;
#' #   1 out of 5 groups test positively with 
#' #   groups of size 10; and
#' #   2 out of 5 groups test positively with 
#' #   groups of size 50.
#' x1 <- c(0,0,1,2)
#' m1 <- c(1,5,10,50)
#' n1 <- c(5,5,5,5)
#' propCI(x=x1, m=m1, n=n1, pt.method="Gart", 
#'        ci.method="skew-score")
#' propCI(x=x1, m=m1, n=n1, pt.method="Gart", 
#'        ci.method="score")
#' 
#' # Reproducing estimates from Table 1 in
#' #   Hepworth & Biggerstaff (2017):
#' propCI(x=c(1,2), m=c(20,5), n=c(8,8), 
#'        pt.method="Firth", ci.method="lrt")
#' propCI(x=c(7,8), m=c(20,5), n=c(8,8), 
#'        pt.method="Firth", ci.method="lrt")

# Brianna Hitt - 03.06.2020
# Removed the "mir" option from pt.method and ci.method arguments
# Removed the "scale" argument from the function (scale = 1 for pooledBin)
# Rewrote parts of the function to make sure that equal group sizes 
#   can be used with other methods since it is just a special case 
#   of unequal group sizes
# Brianna Hitt - 03.09.2020
# 

propCI <- function(x, m, n, pt.method = "mle", 
                   ci.method, conf.level = 0.95, 
                   alternative = "two.sided", 
                   maxiter = 100, tol = .Machine$double.eps^0.5){
  
  if (length(x) != length(m) | length(x) != length(n) | 
      length(m) != length(n)) {
    stop("Arguments x, m, and n must have exactly the same length.\n")
  }
  
  # methods only available for MLE point estimate
  if (pt.method %in% c("Firth", "Gart", "bc-mle") && 
      ci.method %in% c("CP", "Blaker", "AC", "soc", "exact")){
    stop("Please specify a confidence interval method that can be used with the specified point estimate method.\n")
  }
  
  # methods only available for equal group sizes
  if (ci.method %in% c("CP", "Blaker", "AC", "soc") && length(m) > 1){
    stop("Please specify a confidence interval method that can be used when groups are of unequal sizes. See the documentation for valid options.\n")
  }
  
  # determine which confidence interval function to use
  if (ci.method %in% c("skew-score", "score", "bc-skew-score", "lrt", "Wald")) {
    results <- pooledBin(x = x, m = m, n = n, pt.method = pt.method, 
                         ci.method = ci.method, scale = 1, 
                         alpha = 1 - conf.level, tol = tol)
    final <- list("estimate" = results$p, "conf.int" = c(results$lcl, results$ucl), 
                  "pt.method" = pt.method, "ci.method" = ci.method, 
                  "conf.level" = conf.level, "alternative" = "two.sided", 
                  "x" = x, "m" = m, "n" = n)
  } else if (ci.method == "exact") {
    results <- bgtvs(n = n, s = m, y = x, conf.level = conf.level, 
                     alternative = alternative, maxiter = maxiter)
    final <- list("estimate" = results$estimate, "conf.int" = results$conf.int,
                  "pt.method" = "mle", "ci.method" = ci.method, 
                  "conf.level" = conf.level, "alternative" = alternative, 
                  "x" = x, "m" = m, "n" = n)
  } else if (ci.method %in% c("AC", "Blaker", "CP", "soc")) {
    results <- bgtCI(n = n, s = m, y = x, conf.level = conf.level, 
                       method = ci.method, alternative = alternative)
      final <- list("estimate" = results$estimate, "conf.int" = results$conf.int, 
                    "pt.method" = "mle", "ci.method" = ci.method, 
                    "conf.level" = conf.level, "alternative" = alternative, 
                    "x" = x, "m" = m, "n" = n)
  }
  
  class(final) <- "propCI"
  final
}




##################################################################
# print.propCI() function                                        #
##################################################################

# This function was originally written as print.bgtCI and 
#   print.poolbin for binGroup
# Minor modifications were made to combine the two print functions

#' @title Print method for objects of class "propCI"
#' 
#' @description Print method for objects of class "propCI" 
#' created by the \code{\link{propCI}} function.
#' 
#' @param x An object of class "propCI" (\code{\link{propCI}}).
#' @param ... Additional arguments to be passed to \code{print}.
#' 
#' @return A print out of the point estimate and confidence interval 
#' found with \code{\link{propCI}}.
#' 
#' @author This function is a combination of \code{print.poolbindiff} and 
#' \code{print.bgt}, written by Brad Biggerstaff for the \code{binGroup} 
#' package. Minor modifications were made for inclusion of the function in 
#' the \code{binGroup2} package.

# Brianna Hitt - 03.06.2020
# Removed the "scale" argument from the function (scale = 1)

"print.propCI" <- function(x, ...){
  scale <- 1
  
  args <- list(...)
  if (is.null(args$digits)) {digits <- 4}
  else{digits <- args$digits}
  
  switch(x$pt.method, 
         "Firth" = PtEstName <- "Firth's Correction", 
         "Gart" = PtEstName <- "Gart's Correction", 
         "bc-mle" = PtEstName <- "Gart's Correction", 
         "mle" = PtEstName <- "Maximum Likelihood", 
         "mir" = PtEstName <- "Minimum Infection Rate")
  
  switch(x$ci.method, 
         "CP" = CIEstName <- "Clopper-Pearson", 
         "Blaker" = CIEstName <- "Blaker", 
         "AC" = CIEstName <- "Agresti-Coull", 
         "score" = CIEstName <- "Wilson score", 
         "soc" = CIEstName <- "Second-Order Corrected",
         "Wald" = CIEstName <- "Wald", 
         "exact" = CIEstName <- "Exact", 
         "skew-score" = CIEstName <- "Skew-Corrected score (Gart)", 
         "bc-skew-score" = CIEstName <- "Bias- & Skew-Corrected score (Gart)", 
         "lrt" = CIEstName <- "Likelihood Ratio Test Inversion", 
         "mir" = CIEstName <- "Minimum Infection Rate")
  
  cat("\n")
  cat(x$conf.level*100, "percent", CIEstName,
      "confidence interval:\n", sep = " ")
  # cat(scale*x$conf.level*100, "percent confidence interval:\n", sep = " ")
  cat(" [", paste(signif(scale*x$conf.int, digits), collapse = ", "), 
      "]\n", sep = " ")
  cat("Point estimate (", PtEstName, "): ", 
      signif(scale*x$estimate, digits), "\n", sep = "")
  # cat("Point estimate:", signif(scale*x$estimate, digits), "\n", sep = " ")
  # cat("Scale:", scale, sep = " ")
  
  # cat(paste("Point estimator:", PtEstName, "\n"))
  # cat(paste("CI method:", CIEstName, "\n"))
  # cat(paste("Number of individuals:", sum(x$n * x$m), "\n"))
  # cat(paste("Number of groups:", sum(x$n), "\n"))
  # cat(paste("Number of positive groups:", sum(x$x), "\n"))
  
  invisible(x)
}

#

