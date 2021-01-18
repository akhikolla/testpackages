
simgc <- function(locs, sim.n = 1,  marginal, corr, longlat = FALSE){

  if(!is.matrix(locs) & !is.data.frame(locs))
    stop(" Input 'locs' must be a matrix or data frame!")
  ans <- list()
  if(longlat == FALSE){
    D <- as.matrix(dist(locs, method = "euclidean", diag = TRUE, upper = TRUE))
  }else if(longlat == TRUE){
    if (requireNamespace("sp", quietly = TRUE)) {
      D <- sp::spDists(x = as.matrix(locs), longlat = TRUE)
    }else{
      stop("Please install {sp} first!")
    }
  }
  R <- corr$S(D)
  if (requireNamespace("MASS", quietly = TRUE)) {
    sim.grf <- MASS::mvrnorm(sim.n+1, mu = rep(0, nrow(locs)), Sigma = R)
  }else{
    stop("Please install {MASS} first!")
  }
  ans$data <- t(apply(sim.grf, 1, function (x) marginal$q(pnorm(x))))[-(sim.n+1),]
  ans$locs <- locs
  class(ans) <- c("simgc")
  return(ans)
}



#### R Callable Functions
FHUBZIP_R <- function(m1, m2, od1, od2){
  .Call('FHUBZIP', PACKAGE = 'gcKrig', m1, m2, od1, od2)
}

FHUBNB2_R <- function(m1, m2, od1, od2){
  .Call('FHUBNB2', PACKAGE = 'gcKrig', m1, m2, od1, od2)
}

FHUBbinom_R <- function(m1, m2, n1, n2){
  .Call('FHUBbinom', PACKAGE = 'gcKrig', m1, m2, n1, n2)
}


FHUBZIPNB2_R <- function(zipmu, nbmu, zipod, nbod){
  .Call('FHUBZIPNB2', PACKAGE = 'gcKrig', zipmu, nbmu, zipod, nbod)
}


FHUBZIPbinomial_R <- function(zipmu, bmu, zipod, bn){
  .Call('FHUBZIPbinomial', PACKAGE = 'gcKrig', zipmu, bmu, zipod, bn)
}

FHUBNB2binomial_R <- function(nbmu, bmu, nbod, bn){
  .Call('FHUBNB2binomial', PACKAGE = 'gcKrig', nbmu, bmu, nbod, bn)
}



FHUBdiscrete <- function(marg1, marg2, mu1, mu2, od1 = 0, od2 = 0,
                                 binomial.size1 = 1, binomial.size2 = 1)
{
  if (is.function(marg1) | is.function(marg2)){
    stop("Use method 'integral' or 'mc' !")
  }else{
    if (!marg1 %in% c("poisson", "zip", "nb", "binomial"))
      stop("'marg1' must be one of the following: 'poisson', 'zip', 'nb' or 'binomial'. ")
    if (!marg2 %in% c("poisson", "zip", "nb", "binomial"))
      stop("'marg2' must be one of the following: 'poisson', 'zip', 'nb' or 'binomial'. ")

    if(marg1 == 'poisson') od1 = 0; if(marg2 == 'poisson') od2 = 0

    if( (marg1 == 'zip' & marg2 == 'zip') | (marg1 == 'poisson' & marg2 == 'poisson')
        | (marg1 == 'poisson' & marg2 == 'zip') | (marg1 == 'zip' & marg2 == 'poisson') )
      corr <- try(FHUBZIP_R(m1 = mu1, m2 = mu2, od1 = od1, od2 = od2), silent = T)

    if(marg1 == 'nb' & marg2 == 'nb')
      corr <- try(FHUBNB2_R(m1 = mu1, m2 = mu2, od1 = od1, od2 = od2), silent = T)

    if(marg1 == 'binomial' & marg2 == 'binomial')
      corr <- try(FHUBbinom_R(m1 = mu1, m2 = mu2, n1 = binomial.size1,
                              n2 = binomial.size2), silent = T)

    if( (marg1 == 'nb' & marg2 == 'zip') | (marg1 == 'nb' & marg2 == 'poisson') )
      corr <- try(FHUBZIPNB2_R(zipmu = mu2, nbmu = mu1, zipod = od2, nbod = od1), silent = T)
    if( (marg2 == 'nb' & marg1 == 'zip') | (marg2 == 'nb' & marg1 == 'poisson') )
      corr <- try(FHUBZIPNB2_R(zipmu = mu1, nbmu = mu2, zipod = od1, nbod = od2), silent = T)

    if( (marg1 == 'binomial' & marg2 == 'zip') | (marg1 == 'binomial' & marg2 == 'poisson') )
      corr <- try(FHUBZIPbinomial_R(zipmu = mu2, bmu = mu1, zipod = od2, bn = binomial.size1), silent = T)
    if( (marg2 == 'binomial' & marg1 == 'zip') | (marg2 == 'binomial' & marg1 == 'poisson') )
      corr <- try(FHUBZIPbinomial_R(zipmu = mu1, bmu = mu2, zipod = od1, bn = binomial.size2), silent = T)

    if( (marg1 == 'nb' & marg2 == 'binomial'))
      corr <- try(FHUBNB2binomial_R(nbmu = mu1, bmu = mu2, nbod = od1,
                                    bn = binomial.size2), silent = T)
    if( (marg2 == 'nb' & marg1 == 'binomial'))
      corr <- try(FHUBNB2binomial_R(nbmu = mu2, bmu = mu1, nbod = od2,
                                    bn = binomial.size1), silent = T)

    if(corr > 99)
     stop("Marginal mean too large or n too large in binomial:
           computation of the summation is time consuming! Try to use other method.")

    return(corr)
  }
}



corrTG <- function(marg1, marg2, corrGauss = 0.5, method = "integral", nrep = 1000){

 if(!inherits(marg1, "marginal.gc"))  stop("'marg1' must be a function of the class marginal.gc")
 if(!inherits(marg2, "marginal.gc")) stop("'marg2' must be a function of the class marginal.gc")

  atmp <- rep(0,50)
  if(method == "integral"){
    for(u in 1:50){
      atmp[u] = marg1$int.marg(order = u)[[1]]*marg2$int.marg(order = u)[[1]]*(corrGauss)^(u)/factorial(u)
      if (abs(atmp[u]) < .Machine$double.eps^0.5*sqrt(marg1$margvar*marg2$margvar) )
        break
    }
    corr <- sum(atmp)/sqrt(marg1$margvar*marg2$margvar)
  }

  if(method == "mc" & corrGauss == 1){
    corr <- cor(sort(marg1$rsample(nrep = nrep)), sort(marg2$rsample(nrep = nrep)))
  }

  if(method == "mc" & corrGauss < 1){
    if (requireNamespace("MASS", quietly = TRUE)) {
            simnorm = pnorm(MASS::mvrnorm(nrep, mu = c(0,0),
                                          Sigma = matrix(c(1,corrGauss,corrGauss,1),2,2)))
    }else{
      stop("Please install {MASS} first!")
    }
  simmarg1 <- marg1$q(p = simnorm[,1])
  simmarg2 <- marg2$q(p = simnorm[,2])
  corr <- cor(simmarg1, simmarg2)
  }
  corr <- ifelse(corr>1, 1, corr)
  return(corr)
}
