#############################################################################
## Permutation test
#############################################################################
#' Title
#'
#' @param x
#' @param y
#' @param M
#'
#' @return
#' @export
#'
#' @examples
CriticalValue <- function(x, y, sigma, alpha, level, M = 1000, method){
  if(missing(alpha))
    alpha=1
  if(missing(sigma))
    sigma=1  
  if(missing(level))
    level=0.05  
  n <- nrow(as.matrix(x))
  if (length(y)!=n)
    stop( "x and y must be the same size")
  Rho <-rep(0,M)
  if(method=='dCor'| method=='dCov'|method=='gCov' | method=='gCor' | method=='RcppgCov' | method=='RcppgCor'){
    for (i in 1:M){
      y<-y[sample.int(n, replace = FALSE)]
      f <-match.fun(method)
      Rho[i] <- f(x,y, alpha)
    }
    CriticalValue = Rho[order(Rho)][round((1-level)*M)]
  } else if (method=='KdCor'| method=='KdCov'|method=='KgCov' | method=='KgCor' | method=='RcppKgCov' | method=='RcppKgCor'){
    for (i in 1:M){
      y<-y[sample.int(n, replace = FALSE)]
      f <-match.fun(method)
      Rho[i] <- f(x,y, sigma)
    }
    CriticalValue = Rho[order(Rho)][round((1-level)*M)]
  } else {
    stop( "Unrecognized method!")
  }
  return(CriticalValue)
}

#' Title
#'
#' @param x
#' @param y
#' @param M
#' @param sigma
#'
#' @return
#' @export
#'
#' @examples
PermutationTest <- function(x, y, method, sigma, alpha, M = 200, level = 0.05){
  if(missing(alpha))
    alpha=1
  if(missing(sigma))
    sigma=1  
  if(missing(level))
    level=0.05
  n <- nrow(as.matrix(x))
  if (length(y)!=n)
    stop( "x and y must be the same size")
  power = 0
  Rho <-rep(0,M)
  if(method=='dCor'| method=='dCov'|method=='gCov' | method=='gCor' | method=='RcppgCov' | method=='RcppgCor'){
    f <-match.fun(method)
    testStatistics <- f(x,y, alpha)
    for (i in 1:M){
      y<-y[sample.int(n, replace = FALSE)]
      Rho[i] <- f(x,y, alpha)
      if(testStatistics > Rho[i])
        power = power + 1
    }
  } else if (method=='KdCor'| method=='KdCov'|method=='KgCov' | method=='KgCor' | method=='RcppKgCov' | method=='RcppKgCor'){
    f <-match.fun(method)
    testStatistics <- f(x,y, sigma)
    for (i in 1:M){
      y<-y[sample.int(n, replace = FALSE)]
      Rho[i] <- f(x,y, sigma)
      if(testStatistics > Rho[i])
        power = power + 1
    }
    
  } else {
    stop( "Unrecognized method!")
  }
  Pvalue = 1-power/M
  CriticalValue = Rho[order(Rho)][round((1-level)*M)]
  decision = (Pvalue>level)
  if (decision == 0)
    message("Permutation test of method has the critical value of ", CriticalValue, ", the P-value of ", Pvalue, " at significant level ", level, " and we reject the null hypothesis." )
  else
    message("Permutation test of method has the critical value of ", CriticalValue, ", the P-value of ", Pvalue, " at significant level ", level, " and we can not reject the null hypothesis." )    
  return (list(CriticalValue=CriticalValue, Pvalue=Pvalue, method=method, decision=decision))
}
