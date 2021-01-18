# R function for Kcov(x,F(y))
# x, y are same length vectors
# y is order of y to save time
# Kgmd(x,x) computes gini mean difference of x
# KgCor(x,y) = Kgmd(x,y)/Kgmd(x,x)
#' Title
#'
#' Title
#'
#' @param x
#' @param y
#' @param sigma
#'
#' @return
#' @export
#'
#' @examples
Kgmd = function(x, sigma){
  if(is.null(dim(x))) {
    n = length(x)
    Kcovgxy =sum(sqrt(2-2*exp(-as.matrix(dist(as.matrix(x[order(x)])))/sigma)))/(n*(n-1))
    return(Kcovgxy)
  } else{
    x<-as.matrix(x)
    n <- nrow(x)
    return(ifelse(n==1,sum(sqrt(2-2*exp(-as.matrix(dist(as.matrix(x)))/sigma))), sum(sqrt(2-2*exp(-as.matrix(dist(as.matrix(x)))/sigma)))/(n*(n-1))))
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
RcppKGmd = function(x, sigma){
  if(is.null(dim(x))) {
    return(rcpp_Kcovg(x, order(x), sigma))
  } else{
    return(Rcpp_KCovg(x, order(x[,1]), sigma))
  }
}

