#
#    Copyright (C) 2020 David Preinerstorfer
#    david.preinerstorfer@ulb.ac.be
#
#    This file is a part of wbsd.
#
#    wbsd is free software; you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation; either version 2 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
#    GNU General Public License for more details. A copy may be obtained at
#    http://www.r-project.org/Licenses/

###########################################################################
###########################################################################
#
# Compute a bootstrap p value for each column vector of a matrix y.
#
#
###########################################################################
###########################################################################

boot.pval <- function(
#INPUT
y,                            #matrix (n rows) of observations (lhs)
X,                            #design matrix (n times k, rank k)
R,                            #restriction matrix (q times k, rank q)
r,                            #right-hand side in restriction (q-vector)
hcmethod,                     #-1:4; chi2-Test (no df adjustmen), HC0-HC4
restr.cov,                    #Covariance Matrix estimator based on restriction (TRUE or FALSE)
wilddist,                     #Rademacher or Mammen errors in bootstrap samples (only used if comp.meth = "exact")
wildmult,                     #Multiplicator used in bootstrap errors, 0:4 HC0-HC4 weights
wildmult.restr,               #bootstrap multipliers based on restricted hat matrix (TRUE or FALSE)
boot.res.restr,               #bootstrap sample residuals based on restricted res (TRUE or FALSE)
boot.center.restr,            #bootstrap sample centered at restricted res (TRUE or FALSE)
tol = 1e-07,                  #tolerance parameter used in checking invertibility of the covariance matrix in F stat
comp.meth = "exact",          #computation of theta ``exact'' or ``approximation''
Boot.supp = NULL,             #Support bootstrap sample (matrix, n rows, columns = number of bootstrap samples) 
                              #only used if comp.meth = ``approximation''; 
                              #Boot.supp is subsequently further multiplied by wildmult weights
                              #can be generated via auxiliary function: XI.generator
checks = TRUE,                #run input checks (TRUE or FALSE)
cores = 1					            #maximal number of cores used in the computations
)
{

###########################################################################
# Run checks if checks == TRUE
###########################################################################

if(checks == TRUE){
check.y(y, X)
check.X.R(X, R)
check.r(r, R)
check.hcmethod(hcmethod)
check.restr.cov(restr.cov)
check.wilddist(wilddist)
check.wildmult(wildmult)
check.wildmult.restr(wildmult.restr)
check.boot.res.restr(boot.res.restr)
check.tol(tol)
check.comp.meth(comp.meth)
check.Boot.supp(Boot.supp, comp.meth, dim(X)[1])
check.checks(checks)
check.cores(cores)
}

###########################################################################
# Elementary quantities
###########################################################################

n <- dim(X)[1]                        #sample size
k <- dim(X)[2]                        #number of regressors
q <- length(r)                        #number of restrictions
e0 <- matrix(0, nrow = n, ncol = 1)   #vector of zeros of length n
qrX <- qr(X)                          #qr decomposition of X
qrM0lin <- qr(M0lin(X, R))            #qr decomp of basis of M0lin = M0-mu0
Bfac <- Bfactor.matrix(qrX, n, R)     #R(X'X)^{-1}X'
Bfac2 <- Bfactor.matrix2(qrX, R)      #(R(X'X)^{-1}R')^{-1}

###########################################################################
# Prepare input RF for function res.OLSRF (restricted OLS est + resid) 
# [only if needed]
# RF = (X'X)^(-1)R'(R(X'X)^(-1)R')^(-1) = (X'X)^(-1)R' %*% Bfac2
# in case q = k, RF = R^{-1}
###########################################################################

  if(boot.res.restr == TRUE | boot.center.restr == TRUE | restr.cov == TRUE){
   if(q < k){
   factor.tmp2 <- tcrossprod( backsolve(qr.R(qrX), diag(k)) )
   RF <- factor.tmp2%*%t(R)%*%Bfac2
   } else {
   RF <- solve(R)
   }
  } else {
   RF <- NULL
  }
  
###########################################################################
# Generate an element of M0: mu0
###########################################################################

  if(max(abs(r)) > 0){
    mu0 <- X%*%qr.solve(R, r) #assign element of M0 
  } else {
    mu0 <- e0 # if r = 0, then zero vector is in M0lin
  }

###########################################################################
#Prepare wild bootstrap multiplier weights
###########################################################################

  #restricted
  if(wildmult.restr == TRUE){
    bootmult <- c(sqrt(wvec(k-q, n, qrM0lin, wildmult)))
  }
  #unrestricted  
  if(wildmult.restr == FALSE){
    bootmult <- c(sqrt(wvec(k, n, qrX, wildmult)))
  }

###########################################################################
# Generate the support set of the bootstrap sample distribution in case
# of exact computation; compute the probabilities, under which
# each support point is attained under a Rademacher or Mammen measure.
###########################################################################

  if(comp.meth == "exact"){
  
  if(n >= 20){  
  warning("Sample size n might be too large for choosing comp.meth = exact.")
  }

  #In case wilddist = Rademacher:
  #Collect the set of all (2^n) elemts of \{-1, 1\}^n in an (2^n x n) matrix 
  
  #In case wilddist = Mammen:
  #Collect the set of all (2^n) elemts of \{-(sqrt(5)-1)/2, (sqrt(5)+1)/2\}^n 
  #in an (2^n x n) matrix  
  
  if( wilddist == "Rademacher" ){
  XI <- as.matrix(expand.grid(replicate(n, c(-1,1), simplify = FALSE)))
  } else {
  XI <- as.matrix(expand.grid(replicate(n, c(-(sqrt(5)-1)/2, (sqrt(5)+1)/2), 
  simplify = FALSE)))  
  }
  
  #Multiply the vectors generated with the bootstrap multiplicators bootmult

  XI <- XI%*%diag(bootmult)

  #Compute the probability of each row in XI under an n-fold product of 
  #Rademacher or Mammen distributions

  #the probabilities of all coordinates in XI under independent
  #Rademacher distributions
  if(wilddist == "rademacher"){
  probs <- expand.grid(replicate(n, c(.5,.5), simplify = FALSE))
  } 

  #the probabilities of all coordinates in XI under independent
  #Mammen distributions
  if(wilddist == "mammen"){
  probs <- expand.grid(replicate(n, c((sqrt(5)+1)/(2*sqrt(5)), 
                      (sqrt(5)-1)/(2*sqrt(5))), simplify = FALSE))
  }

  #combine to get the probability of each row of XI under Rademacher or
  #Mammen distributions
  probs <- apply(probs,1,prod)
  }

###########################################################################
# Generate the support set of the bootstrap sample distribution in case
# of approximate computation; i.e., transpose input Boot.supp, and multiply by
# diag(bootmult)
###########################################################################

  if(comp.meth == "approximation"){
  XI <- t(Boot.supp)%*%diag(bootmult)
  }

###########################################################################
# Computation of the bootstrap p values
###########################################################################

#create storage space for p-values and compute test statistic on data y

  pvals <- rep(NA, length = dim(y)[2])
  Tvals.comp <- F.wrap(y, R, r, r, X, n, k, q, qrX, 
                    hcmethod, cores, restr.cov, 
                    Bfac, Bfac2, qrM0lin, RF, tol)

  for(i in 1:dim(y)[2]){

  #generate the bootstrap sample
  
  Ystar <- boot.sample(y[,i, drop = FALSE], t(XI), qrX, 
              R, r, boot.res.restr, boot.center.restr, RF)
              
    if(boot.center.restr == FALSE){
      rstar <- R%*%qr.coef(qrX, y[,i])
    } else {
      rstar <- r
    }
    
    Tvalsboot <- F.wrap(Ystar, R, r, rstar, X, n, k, q, qrX, 
                    hcmethod, cores, restr.cov, 
                    Bfac, Bfac2, qrM0lin, RF, tol)
                    
    select.support <- ( (Tvalsboot >= Tvals.comp[i]) | (Tvalsboot == -1) )
    if(comp.meth == "approximation"){
      pvals[i] <- mean(select.support)
    } else {
      pvals[i] <- sum((probs[select.support]))
    }
    
    }
return(list("p" = pvals))
}