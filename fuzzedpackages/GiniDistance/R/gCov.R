#' Title
#'
#' @param x
#' @param y
#'
#' @return
#' @export
#'
#' @examples
gCov = function(x,y, alpha){
  if(missing(alpha))
    alpha <-1
  if (alpha==1){
    if(is.null(dim(x))){
      n <- nrow(as.matrix(x))
      m <- ncol(as.matrix(x))
      if (length(y)!=n)
        stop( "x and y must be the same size")
      z<-cbind(x,y)
      z<-z[ order(z[,m+1], z[,1]),]
      delta<-rep(0,n)
      max <- as.numeric(z[1,m+1])
      group <-1
      n1<-1
      n2<-1
      for(i in 2:n){
        if (as.numeric(z[i,2]) == max){
          n2=i
        } else {
          delta[group]=gmd(z[n1:n2,1])*(n2-n1+1)/n
          n1<-i
          max<-as.numeric(z[n2+1,2])
          group=group+1
        }
        if(i==n){
          delta[group]=gmd(z[n1:n2,1])*(n2-n1+1)/n
        }
      }
      gCovxy = gmd(z[,1])-sum(delta)
      return(abs(gCovxy))
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
          delta[group]=gmd(z[n1:n2,1:m])*(n2-n1+1)/n
          n1<-i
          max<-as.numeric(z[n2+1,m+1])
          group=group+1
        }
        if(i==n){
          delta[group]=gmd(z[n1:n2,1:m])*(n2-n1+1)/n
          
        }
      }
      gCovxy = abs(gmd(z[,1:m])-sum(delta))
      return(gCovxy)
    }
  } else if (0<alpha & 2>=alpha){
    x<- as.matrix(x)
    n <- nrow(as.matrix(x))
    m <- ncol(as.matrix(x))
    if (length(y)!=n)
      stop( "x and y must be the same size")
    z<-cbind(x,y)
    z<-z[order(z[,m+1], z[,1]),]
    
    u<-as.matrix(dist(z[,1:m]))^alpha
    
    delta<-rep(0,n)
    
    max <- as.numeric(z[1,m+1])
    group <-1
    n1<-1
    n2<-1
    
    for(i in 2:n){
      if (as.numeric(z[i,m+1]) == max){
        n2=i
      } else {
        delta[group]=sum(u[n1:n2,n1:n2])/((n2-n1)*n)
        n1<-i
        max<-as.numeric(z[n2+1,m+1])
        group=group+1
      }
      if(i==n){
        delta[group]=sum(u[n1:n2,n1:n2])/((n2-n1)*n)
        
      }
    }
    gCovxy = sum(u)/((n-1)*n)-sum(delta)
    return(abs(gCovxy))
  } else{
    stop( "alpha must be the a number between 0 and 2")
  }
  
}


RcppgCov = function(x,y,alpha){
  if(missing(alpha))
    alpha <-1
  if (alpha==1){
    if(is.null(dim(x))){
      n <-length(x)
      m=1
    } else {
      x<- as.matrix(x)
      n <- nrow(x)
      m <- ncol(x)
    }
    if (length(y)!=n)
      stop( "x and y must be the same size")
    z<-cbind(x,y)
    z<-z[ order(z[,m+1], z[,1]),]
    return(abs(Rcpp_gCov(z)))
  } else if (0<alpha & 2>=alpha){
    x<-as.matrix(x)
    n <- nrow(x)
    m <- ncol(x)
    if (length(y)!=n)
      stop( "x and y must be the same size")
    z<-cbind(x,y)
    z<-z[ order(z[,m+1], z[,1]),]
    return(ifelse(n==1, 0, Rcpp_gCov_Alpha(z, alpha)))
  } else{
    stop( "alpha must be the a number between 0 and 2")
  }
}

