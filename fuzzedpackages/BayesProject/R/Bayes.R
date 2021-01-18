# *********************************************************************************************************************
# Rcpp wrapper functions which all return a list with $cusum and $proj and $cpt
# exception: wildBinSeg only returns a changepoint vector
# all functions rescale variance like Inspect so that threshold can be used from compute.threshold()
# *********************************************************************************************************************
#' @useDynLib BayesProject



# not exported
rescale.variance <- function(x) {
  p <- nrow(x)
  for (j in 1:p) {
    scale <- mad(diff(x[j, ]))/sqrt(2)
    x[j, ] <- x[j, ]/scale
  }
  return(x)
}

# not exported
single.change <- function(n, p, k, z, vartheta, sigma = 1, shape = 3, noise = 0) {
  mu <- matrix(0, p, n)
  if (shape == 3) {
    theta = (1:k)^(-1/2)
    theta = theta/norm(theta,type="2") * sqrt(k) * vartheta
  }
  if (noise != -1) {
    mu[1:k, (z + 1):n] <- theta
  }
  if (noise <= 0) {
    W = matrix(rnorm(p * n), p, n) * sigma
  }
  x <- mu + W
  return(list(x=x,mu=mu))
}

# not exported
# function to calibrate the changepoint detection threshold, requires function "cusumFct" to compute cusum(x) with single parameter for the dataset x
compute.threshold <- function(n, p, cusumFct, nrep=100) {
  cusum.stats <- sapply(1:nrep, function(i) { max(cusumFct(single.change(n, p, 1, n - 1, 0)$x)) })
  max(cusum.stats)
}

# not exported
cusum.spacing <- function(n,maxTau,alpha=0.6) {
  # points t, truth tau
  # tau <= t
  t <- maxTau
  v <- t
  for(tau in maxTau:1) {
    if(tau/t<alpha) {
      t <- tau
      v <- c(t,v)
    }
  }
  # t <= tau
  t <- 1
  w <- t
  for(tau in 1:maxTau) {
    if((n-tau)/(n-t)<alpha) {
      t <- tau
      w <- c(w,t)
    }
  }
  return(sort(union(v,w)))
}

# not exported
vector.cusum.transform <- function(x) {
  n <- length(x)
  leftsums <- cumsum(x)
  rightsums <- leftsums[n] - leftsums
  t <- 1:(n - 1)
  return( (rightsums[t]/(n - t) - leftsums[t]/t) * sqrt(t*(n - t)/n) )
}



#' Cpp implementation of the Bayesian projection algorithm to detect single multivariate changepoints.
#' 
#' Detects one multivariate changepoint in a dataset using the fast projection direction algorithm of Hahn et al. (2019).
#' Solely required is the dataset as first parameter.
#' The testing threshold ("threshold"), the number of timepoints to calculate a projection ("nTimePoints") and the regularisation parameter ("K") are chosen automatically.
#' 
#' @param x A \eqn{p \times n} matrix representing \eqn{p} data series having \eqn{n} observations each.
#' @param threshold The testing threshold to detect the single changepoint. If missing, parameter will be calibrated automatically.
#' @param nTimePoints The number of equidistant timepoints at which the projection direction is calculated. If no value (NULL) is given, timepoints are chosen automatically.
#' @param K The regularisation parameter for the Bayesian projection direction. Default is \eqn{1/\sqrt(2)}.
#' @param rescale.var A boolean flag to indicate if the variance should be rescaled before detecting a changepoint. Default is TRUE.
#' 
#' @importFrom Rdpack reprompt
#' @references Hahn, G., Fearnhead, P., Eckley, I.A. (2020). Fast computation of a projection direction for multivariate changepoint detection. Stat Comput.
#' 
#' @examples
#' library(BayesProject)
#' data(testdata)
#' res <- bayes(testdata,nTimePoints=100)
#' print(res$cpt)
#' 
#' @export
bayes <- function(x, threshold, nTimePoints=NULL, K=1/sqrt(2), rescale.var=TRUE) {
  if(rescale.var) x <- rescale.variance(x)
  p <- nrow(x)
  n <- ncol(x)
  if(is.null(nTimePoints)) {
    timePoints <- cusum.spacing(n,n-1)
  }
  else {
    timePoints <- floor(seq(1,n-1,length.out=nTimePoints))
  }
  if(missing(threshold)) threshold <- compute.threshold(n=n,p=p,cusumFct=function(d) bayes_cpt(d,timePoints,K)[-(1:p)])
  
  v <- bayes_cpt(x, timePoints, K)
  proj <- v[1:p]
  # cusum <- v[-(1:p)]
  fullcusum <- abs(vector.cusum.transform(proj %*% x))
  return(list( cusum=fullcusum, proj=proj, cpt=ifelse(max(fullcusum)>threshold,which.max(fullcusum)[1],NA) ))
}



#' Cpp implementation of sum-cusum and max-cusum for single changepoint detection.
#' 
#' Detects one multivariate changepoint in a dataset using the sum-cusum or max-cusum technique.
#' Solely required is the dataset as first parameter.
#' The testing threshold ("threshold") is chosen automatically if missing.
#' The parameter "sum_cusum" (default TRUE) indicates if sum-cusum or max-cusum is used.
#' 
#' @param x A \eqn{p \times n} matrix representing \eqn{p} data series having \eqn{n} observations each.
#' @param threshold The testing threshold to detect the single changepoint. If missing, parameter will be calibrated automatically.
#' @param sum_cusum A boolean flag to indicate if sum cusum (sum_cusum=T) or max cusum (sum_cusum=F) is used. Default is TRUE.
#' @param rescale.var A boolean flag to indicate if the variance should be rescaled before detecting a changepoint. Default is TRUE.
#' 
#' @importFrom Rdpack reprompt
#' @references Hahn, G., Fearnhead, P., Eckley, I.A. (2020). Fast computation of a projection direction for multivariate changepoint detection. Stat Comput.
#' 
#' @examples
#' library(BayesProject)
#' data(testdata)
#' resSumCusum <- summaxcusum(testdata,sum_cusum=TRUE)
#' print(resSumCusum$cpt)
#' resMaxCusum <- summaxcusum(testdata,sum_cusum=FALSE)
#' print(resMaxCusum$cpt)
#' 
#' @export
summaxcusum <- function(x, threshold, sum_cusum=TRUE, rescale.var=TRUE) {
  if(rescale.var) x <- rescale.variance(x)
  p <- nrow(x)
  n <- ncol(x)
  if(missing(threshold)) threshold <- compute.threshold(n=n,p=p,cusumFct=function(d) sum_max_cusum(d,sum_cusum)[-(1:p)])
  
  v <- sum_max_cusum(x, sum_cusum)
  proj <- v[1:p]
  cusum <- v[-(1:p)]
  return(list( cusum=cusum, proj=proj, cpt=ifelse(max(cusum)>threshold,which.max(cusum)[1],NA) ))
}



#' Wild Binary Segmentation Wrapper for the functions "bayes" and "summaxcusum".
#' 
#' Detects multivariate changepoints in a dataset using the Wild Binary Segmentation framework of Fryzlewicz (2014).
#' The dataset is supplied as the first parameter.
#' The second parameter is a calibrated function "cusumFct(x)" which takes a multivariate data matrix \eqn{x} and returns a cusum vector for it.
#' The threshold is supplied with parameter "threshold", parameter "m" specifies the number of random WBS intervals on each recursion level,
#' and "minwindow" is the minimal window size up to which the dataset is further divided recursively to find more changepoints.
#' 
#' @param x A \eqn{p \times n} matrix representing \eqn{p} data series having \eqn{n} observations each.
#' @param cusumFct A calibrated function which returns a cusum vector for a dataset supplied as its single input parameter. Note that rescaling of the variance should be deactivated inside "cusumFct". When using the function "bayes" as in the example below, it is advised to set the threshold to e.g. zero in order to deactivate the time-consuming (and unnecessary) threshold computation inside "bayes".
#' @param threshold The testing threshold to detect the single changepoint. The threshold must be specified.
#' @param m The number of random WBS intervals on each recursion level.
#' @param minwindow The minimal window size up to which the dataset is further divided recursively to find more changepoints.
#' @param rescale.var A boolean flag to indicate if the variance should be rescaled before detecting a changepoint. Default is TRUE.
#' 
#' @importFrom Rdpack reprompt
#' @references Fryzlewicz, P. (2014). Wild binary segmentation for multiple change-point detection. Ann Statist, 42(6):2243--2281.
#' 
#' @examples
#' library(BayesProject)
#' data(testdata)
#' bayes_cusum <- function(x) bayes(x,threshold=0,rescale.var=FALSE)$cusum
#' res <- wildBinSeg(testdata, cusumFct=bayes_cusum, threshold=1)
#' print(res)
#' 
#' @export
wildBinSeg <- function(x, cusumFct, threshold, m=100, minwindow=10, rescale.var=TRUE) {
  if(rescale.var) x <- rescale.variance(x)
  if(!is.matrix(x)) x <- matrix(x,nrow=1)
  n <- ncol(x)
  
  wildbinseg <- function(s,e) {
    if(s>e-minwindow) return(NULL)
    if(s<(e-minwindow)) M <- sample(s:(e-minwindow), size=m, replace=T)
    else M <- rep(s,m)
    M2 <- sapply(M,function(i) ifelse((i+minwindow)<e, sample((i+minwindow):e,size=1), e) )
    M <- cbind(M,M2)
    colnames(M) <- NULL
    # compute maximal cusum statistic for each interval
    allcusum <- sapply(1:m, function(i) max(cusumFct(x[,M[i,1]:M[i,2]])) )
    # changepoint found? -- determine wbs interval in which the changepoint was found, and then its exact position
    if(max(allcusum,na.rm=T)>threshold) {
      int <- which.max(allcusum)
      cpt <- M[int,1]-1 + which.max( cusumFct(x[,M[int,1]:M[int,2]]) )
      return(c( wildbinseg(s,cpt), cpt, wildbinseg(cpt+1,e) ))
    }
    else return(NULL)
  }
  wildbinseg(1,n)
}
