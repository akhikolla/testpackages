k1Fun1Cos <- function(x) {
  x <- as.vector(x)
  res <- cos(x)
  attr(res, "der") <- matrix(NA, nrow = length(x), ncol = 2)
  attr(res, "der")[, 1] <- - sin(x)
  attr(res, "der")[, 2] <- - res
  return(res)
}

k1Fun1Exp <- function(x) {
    ## XXX should we keep some attributes of 'x' as do classical function
    ## such as 'exp' or 'dnorm'???
    ##
    ## d <- attributes(x)$dim
    res <- .Call(k1FunExpC, x)
    ## attr(res, "dim") <- d
    res
}

k1Fun1Matern3_2 <- function(x) {
    res <- .Call(k1FunMatern3_2C, x)
    res
}

k1Fun1Matern5_2 <- function(x) {
    res <- .Call(k1FunMatern5_2C, x)
    res
}

##' Plain one-dimensional kernel function.
##'
##' @title Plain One-Dimensional Kernel Function.
##' 
##' @param x A vector of numeric values.
##' 
##' @return A numeric vector with an attribute named \code{"der"}
##' containing the first two derivatives of the function.
##'
##' @note There can be some changes in the future concerning the way
##' of coping with the attributes of \code{'x'}. If \code{x}
##' is a matrix or array, should the returned value and its attribute
##' be an array too. Should dimnames be kept?
##' 
k1Fun1Gauss <- function(x) {
    res <- .Call(k1FunGaussC, x)
    res
}

k1Fun1PowExp <- function(x, alpha = 1.5) {
    res <- .Call(k1FunPowExpC, x, alpha)
    res
}
