# R function for cov(x,F(y))
# x, y are same length vectors
# y is order of y to save time
# covg(x,x) computes gini mean difference of x
# corg(x,y) = covg(x,y)/covg(x,x)
#' Title
#'
#' @param x
#' @param y label
#'
#' @return Gini covariance
#' @export
#'
#' @examples
gmd = function(x, alpha){
  if(missing(alpha)) 
    alpha <- 1
  if(alpha==1){
    if(is.null(dim(x))) {
      n = length(x)
      r = 2*seq(1,n)-n-1
      covgxy =(2*sum(r*x[order(x)])/(n*(n-1)))
      return(abs(covgxy))
    } else{
      x<-as.matrix(x)
      n <- nrow(x)
      x <- as.matrix(dist(x))
      return(ifelse(n==1, 0, abs(sum(x)/(n*(n-1)))))
    }
  } else if (0<alpha & 2>=alpha){
    x<-as.matrix(x)
    n <- nrow(x)
    x <- as.matrix(dist(x))
    return(ifelse(n==1, 0, abs(sum(x^alpha)/(n*(n-1)))))
    
  } else{
    stop( "alpha must be the a number between 0 and 2")
  }
}

# R function for cov(x,F(y))
# x, y are same length vectors
# y is order of y to save time
# covg(x,x) computes gini mean difference of x
# corg(x,y) = covg(x,y)/covg(x,x)
#' Title
#'
#' @param x
#' @param y label
#'
#' @return Parallel Rcpp Gini covariance
#' @export
#'
#' @examples
RcppParallelGmd = function(x){
  if(is.null(dim(x))) {
    return(rcpp_covg(x, order(x)))
  } else{
    return(Rcpp_Parallel_Covg(x, order(x)))
  }
}
# R function for cov(x,F(y))
# x, y are same length vectors
# y is order of y to save time
# covg(x,x) computes gini mean difference of x
# corg(x,y) = covg(x,y)/covg(x,x)
#' Title
#'
#' @param x
#' @param y label
#'
#' @return Rcpp Gini covariance
#' @export
#'
#' @examples
RcppGmd = function(x, alpha){
  if(missing(alpha)) 
    alpha <- 1
  if(alpha==1){
    if(is.null(dim(x))) {
      return(rcpp_covg(x, order(x)))
    } else{
      return(Rcpp_Covg(x,order(x)))
    }
  } else if (0<alpha & 2>=alpha){
    x<-as.matrix(x)
    n <- nrow(x)
    m <- ncol(x)
    if(m==1){
      return(ifelse(n==1, 0, rcpp_covg_alpha(x, order(x), alpha)))
    } else {
      return(ifelse(n==1, 0, Rcpp_Covg_Alpha(x, order(x), alpha)))
    }
  } else{
    stop( "alpha must be the a number between 0 and 2")
  }
  
}
