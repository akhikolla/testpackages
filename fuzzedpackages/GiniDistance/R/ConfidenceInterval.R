#' Title
#'
#' @param x
#' @param y
#'
#' @return
#' @export
#'
#' @examples
ConfidenceInterval <- function(x, y, sigma, alpha, level, method){
  if(missing(alpha))
    alpha=1
  if(missing(sigma))
    sigma=1
  if (missing(level) | level<0 | level>1)
    level = 0.95
  n <- nrow(as.matrix(x))
  m <- ncol(as.matrix(x))
  if (length(y)!=n)
    stop( "x and y must be the same size")
  rhoHat <-rep(0,n)
  if(method=='dCor'| method=='dCov'|method=='gCov' | method=='gCor' | method=='RcppgCov' | method=='RcppgCor'){
    if(m==1){
      for(i in 1:n){
        x1<-x[-i]
        y1<-y[-i]
        f <-match.fun(method)
        rhoHat[i]<-f(x1,y1, alpha)
      }
    }
    else{
      for(i in 1:n){
        x1<-x[-i,]
        y1<-y[-i]
        f <-match.fun(method)
        rhoHat[i]<-f(x1,y1, alpha)
      }
    }
    mean = f(x,y, alpha)
    s = sd(rhoHat)*(n-1)/sqrt(n)
    z = qnorm(0.5+level/2)
    return(c(mean-z*s, mean+z*s))
  } else if (method=='KdCor'| method=='KdCov'|method=='KgCov' | method=='KgCor' | method=='RcppKgCov' | method=='RcppKgCor'){ 
    if(m==1){
      for(i in 1:n){
        x1<-x[-i]
        y1<-y[-i]
        f <-match.fun(method)
        rhoHat[i]<-f(x1,y1, sigma)
      }
    }
    else{
      for(i in 1:n){
        x1<-x[-i,]
        y1<-y[-i]
        f <-match.fun(method)
        rhoHat[i]<-f(x1,y1, sigma)
      }
    }
    mean = f(x,y, alpha)
    s = sd(rhoHat)*(n-1)/sqrt(n)
    z = qnorm(0.5+level/2)
    return(c(mean-z*s, mean+z*s))
  } else{
    stop( "Unrecognized method!")
  }
}

