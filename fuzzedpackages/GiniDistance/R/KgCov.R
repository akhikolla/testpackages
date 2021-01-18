#' Title Kernelized Gini covariance for multivariate X and categorical Y ###
#'
#' @param x
#' @param y
#' @param sigma
#'
#' @return
#' @export
#'
#' @examples
KgCov = function(x,y, sigma){
    if(is.null(dim(x))){
      n <- nrow(as.matrix(x))
      if (length(y)!=n)
        stop( "x and y must be the same size")
      z<-cbind(x,y)
      
      z<-z[ order(z[,2], z[,1]),]
      
      delta<-rep(0,n)
      
      max <- as.numeric(z[1,2])
      group <-1
      n1<-1
      n2<-1
      
      for(i in 2:n){
        if (as.numeric(z[i,2]) == max){
          n2=i
        } else {
          delta[group]=Kgmd(z[n1:n2,1], sigma)*(n2-n1+1)/n
          n1<-i
          max<-as.numeric(z[n2+1,2])
          group=group+1
        }
        if(i==n){
          delta[group]=Kgmd(z[n1:n2,1], sigma)*(n2-n1+1)/n
          
        }
      }
      KgCovxy = Kgmd(x, sigma)-sum(delta)
      return(abs(KgCovxy))
    } else{
      x<- as.matrix(x)
      n <- nrow(as.matrix(x))
      m <- ncol(as.matrix(x))
      if (length(y)!=n)
        stop( "x and y must be the same size")
      z<-cbind(x,y)
      z<-z[order(z[,m+1], z[,1]),]
      
      delta<-rep(0,n)
      
      max <- as.numeric(z[1,m+1])
      group <-1
      n1<-1
      n2<-1
      
      for(i in 2:n){
        if (as.numeric(z[i,m+1]) == max){
          n2=i
        } else {
          delta[group]=Kgmd(z[n1:n2,1:m], sigma)*(n2-n1+1)/n
          n1<-i
          max<-as.numeric(z[n2+1,m+1])
          group=group+1
        }
        if(i==n){
          delta[group]=Kgmd(z[n1:n2,1:m], sigma)*(n2-n1+1)/n
          
        }
      }
      KgCovxy = abs(Kgmd(z[,1:m], sigma)-sum(delta))
      return(KgCovxy)
    }
}


#' Title Kernelized Gini covariance for multivariate X and categorical Y ###
#'
#' @param x
#' @param y
#' @param sigma
#'
#' @return
#' @export
#'
#' @examples
RcppKgCov = function(x,y, sigma){
    if(is.null(dim(x))){
      n <- length(x)
      m =1
    } else{
      x<- as.matrix(x)
      n <- nrow(as.matrix(x))
      m <- ncol(as.matrix(x))
    }
    if (length(y)!=n)
      stop( "x and y must be the same size")
    z<-cbind(x,y)
    z<-z[ order(z[,m+1], z[,1]),]
    return(Rcpp_KgCov(z, sigma))
}
