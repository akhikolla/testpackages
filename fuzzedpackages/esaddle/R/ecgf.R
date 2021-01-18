#####
#' Cumulant generating function estimation
#' @description Calculates the empirical cumulant generating function (CGF) and its derivatives
#'               given a sample of n d-dimentional vectors.
#'
#' @param lambda point at which the empirical CGF is evaluated (d-dimensional vector).
#' @param X an n by d matrix containing the data.
#' @param mix fraction of empirical and normal CGF to use. If \code{mix==1} only the empirical CGF is used, 
#'        if \code{mix==0} only the normal CGF is used.
#' @param grad if \code{grad==0} only the value of the CGF at \code{lambda} is returned, 
#'             if \code{grad==1} also its first derivative wrt \code{lambda} 
#'             and \code{if grad==2} also the second derivarive wrt \code{lambda}.
#' @return A list with entries:
#'         \itemize{
#'         \item{ \code{K} }{the value of the empirical CGF at \code{lambda};}
#'         \item{ \code{dK} }{the value of the gradient empirical CGF wrt lambda at \code{lambda};}
#'         \item{ \code{d2K} }{the value of the hessian of the empirical CGF wrt lambda at \code{lambda}.}
#'         }
#' @details For details on the CGF estimator being used here, see Fasiolo et al. (2016).
#' @author Matteo Fasiolo <matteo.fasiolo@@gmail.com> and Simon N. Wood.
#' @references Fasiolo, M., Wood, S. N., Hartig, F. and Bravington, M. V. (2016). 
#'             An Extended Empirical Saddlepoint Approximation for Intractable Likelihoods. ArXiv http://arxiv.org/abs/1601.01849.
#' @examples 
#' X <- matrix(rnorm(2 * 1e3), 1e3, 2)
#' K <- ecgf(lambda = c(0, 0), X = X, mix = 0.5, grad = 2) 
#' K$K # CGF
#' K$dK # CGF' (gradient)
#' K$d2K # CGF'' (Hessian)
#' @export
#'
ecgf <- function(lambda, X, mix, grad = 0) {
  ## X[i,j] is ith rep of jth variable. Evaluate observed KGF 
  ## and its derivs w.r.t. lambda, without overflow...
  
  out <- .ecgf(lambda = lambda, 
               X = X, 
               kum1 = colMeans(X), 
               kum2 = .robCov(t(X), alpha2 = 4, beta2 = 1.25)$COV, 
               grad = grad, 
               mix = mix )
  
  return( out[c("K", "dK", "d2K")] ) 
  
}

####
# Same as ecgf() but this require also kum1, kum2, lets you choose the mixMethod 
# and returns the also the original estimates of K, K', K'', together with their tilted versions
#
.ecgf <- function(lambda, X, kum1, kum2, mix, grad) {
  ## X[i,j] is ith rep of jth variable. Evaluate observed KGF 
  ## and its derivs w.r.t. lambda, without overflow...
  
  if(!is.vector(lambda)) lambda <- as.vector(lambda) 
  if (!is.matrix(X)) X <- matrix(X, length(X), 1)
  n <- nrow(X)
  d <- ncol(X)
  
  stopifnot(d == length(lambda))
  
  ret <- .Call("ecgfCpp",
               lambda_ = lambda, 
               X_ = X, 
               mix_ = mix, 
               grad_ = grad, 
               kum1_ = kum1, 
               kum2_ = kum2)
    
  return( ret )
  
}




