#' Jackknife Mutual Information
#'
#' This function provides method for dependence test. It uses permutation test to determine
#' the rejection region.
#' @param x n by p sample matrix.
#' @param y n by q sample matrix.
#' @param BN Number of permutations, the default value is 1000.

#' @return the output is a list which contains:
#' \itemize{
#' \item mi: the value of Jackknife Mutual information
#' \item pvalue: the p-value of independence test that based on the permutation of JMI, the value is not provided if BN=0.
#' }
#'
#'
#' @export
#'
#' @references Zeng, X., Xia, Y., & Tong, H. (2018). Jackknife approach to the estimation of mutual information[J]. Proceedings of the National Academy of Sciences, 201715593.
#'
#' @examples
#'  x <- matrix(rnorm(50*3),50,3)
#'  y <- matrix(rnorm(50*2),50,2)
#'  #calculate the Jackknife Mutual information between x and y.
#'  JMI(x,y,0)$mi
#'  #calculate the p-value of independent test between x and y that based on 500 permutations.
#'  JMI(x,y,500)$pvalue
#'
#'
#' @useDynLib JMI
#' @import Rcpp



JMI<-function(x,y,BN=1000){
  x=as.matrix(x);
  y=as.matrix(y);
  if (BN>0)
  {
  res <- mJMICpp(x,y,BN)
  return(res)
  }
  if (BN==0)
  {
    res <- mJMICpp(x,y,BN)$mi
    return(list(mi=res))
  }
}
