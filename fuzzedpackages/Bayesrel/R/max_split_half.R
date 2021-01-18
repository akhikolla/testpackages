# this function uses an exhaustive search
# the latter is much faster when the number of items increases e.g. 10 and up
# source:
# Benton, T. (2015). An empirical assessment of Guttman’s Lambda 4 reliability coefficient.
# In Quantitative Psychology Research (pp. 301-310). Springer, Cham.
#

MaxSplitExhaustive <- function(M){
  #data – matrix of items scores (row=candidates,column=items)
  cov1 <- M
  nite <- ncol(M)
  mat1 <- (bin.combs2(nite)+1)/2
  res1 <- 0
  for (jjz in 1:length(mat1[,1])){
    xal <- mat1[jjz,]
    gutt1 <- 4*(t(xal)%*%cov1%*%(1-xal))/sum(cov1)
    resrand <- gutt1
    if (resrand > res1){
      res1 <- resrand
    }
  }
  return(res1)
}

bin.combs2 <- function (p) {
  retval <- matrix(0, nrow = 2^p, ncol = p)
  for (n in 1:p) {
    retval[, n] <- rep(c(rep(-1, (2^p/2^n)), rep(1, (2^p/2^n))),
                       length = 2^p)
  }
  len <- (nrow(retval)/2)
  comb <- retval[1:len, ]
  comb
}


# finds split reliability coefficients, produces sample of lambda4s,
# default return is max lambda4, but other quantiles can also be extracted
# source:
# Hunt, T. D., & Bentler, P. M. (2015). Quantile lower bounds to reliability based on locally optimal splits.
# Psychometrika, 80(1), 182-195.
quant.lambda4 <- function(x, starts = 1000, quantile = 1){
  l4.vect <- rep(NA, starts)
  #Determines if x is a covariance or data matrix and establishes a covariance matrix for estimation.
  sigma  <- x
  items <- ncol(sigma)
  #Creates an empty matrix for the minimized tvectors
  splitmtrx <- matrix(NA, nrow = items, ncol = starts)
  # creates the row and column vectors of 1s for the lambda4 equation.
  onerow <- rep(1,items)
  onerow <- t(onerow)
  onevector <- t(onerow)
  f <- rep(NA,starts)
  for(y in 1:starts){
    #Random number generator for the t-vectors
    trow <- (round(runif(items, min  = 0, max = 1))-.5)*2
    trow <- t(trow)
    tvector <- t(trow)
    #Creating t vector and row
    tk1 <- (tvector)
    tk1t <- t(tk1)
    tk2 <- (trow)
    tk2t <- t(tk2)
    #Decision rule that determines which split each item should be on. Thus minimizing the numerator.
    sigma0 <- sigma
    diag(sigma0) <- 0
    random.order <- sample(1:items)
    for (o in 1:items){
      oi <- sigma0[,random.order[o]]
      fi <- oi%*%tk1
      if (fi < 0) {tk1[random.order[o],1] <- 1}
      if (fi >= 0) {tk1[random.order[o],1] <- -1}
    }
    t1 <- (1/2)*(tk1+1)
    fk1 <- tk1t%*%sigma0%*%tk1
    t1t <- t(t1)
    t2 <- 1*(1-t1)
    t2t <- t(t2)
    f[y] <- fk1
    splitmtrx[,y] <- t1
    l4.vect[y] <- (4*(t1t%*%sigma%*%t2))/(onerow%*%sigma%*%onevector)
  }
  quants <- quantile(l4.vect, quantile)
  lambda4.quantile <- quants
  return(lambda4.quantile)
}

