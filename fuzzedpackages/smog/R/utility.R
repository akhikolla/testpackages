# Copyright (c) 2018 - 2020 Chong Ma
 
# This file contains the kernel function for the R package smog. 
# The function smog is written for the generalized linear model constraint 
# on specified hierarchical structures by using overlapped group penalty. 
# It is implemented by combining the ISTA and ADMM algorithms, and works 
# for continuous, multimonial and survival data. 


#' smog generic
#'
#' @param x an object of class from ''smog''.
#' @param ... further arguments passed to or from other methods.
#' @keywords internal
#' 
#' @seealso \code{\link{cv.smog}}, \code{\link{smog.default}}, \code{\link{smog.formula}}, 
#'          \code{\link{predict.smog}}, \code{\link{plot.smog}}.
#' @return the cofficients table of the object x of class from ''smog''.
#'         See \code{\link[base]{print}}.
#' 
#' @references \insertRef{ma2019structural}{smog}
#' 
#' @export
smog <- function(x, ...) UseMethod("smog")


#' @rdname smog
#' 
#' @keywords internal
#' 
#' @export
print.smog <- function(x, ...){
  cat("Call:\n")
  print(x$call)
  
  cat("\nCoefficients:\n")
  print(x$coefficients,row.names=FALSE)
}


#' @rdname cv.smog
#' 
#' @keywords internal
#' 
#' @export
print.cv.smog <- function(x, ...){
  cat("Call:\n")
  print(x$call)
  
  cat("\nCoefficients:\n")
  print(x$cvfit)
}


#' @rdname cv.cglasso
#' 
#' @keywords internal
#' 
#' @export
print.cv.cglasso <- function(x, ...){
  cat("Call:\n")
  print(x$call)
  
  cat("\nCoefficients:\n")
  print(x$cvfit)
}





