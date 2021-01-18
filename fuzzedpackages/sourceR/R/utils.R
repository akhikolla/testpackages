#####################################################
# Name: utils.R                                     #
# Author: Chris Jewell <c.jewell@lancaster.ac.uk>   #
# Created: 20161206                                 #
# Copyright: Chris Jewell 2016                      #
# Purpose: Miscellaneous helper functions           #
#####################################################


isFiniteInteger = function(x)
{
  int = x
  mode(x) = 'integer'
  is.finite(x) & isTRUE(all.equal(x, int))
}
isFiniteLogical = function(x)
{
  is.finite(x) & is.logical(x)
}
isFinitePositive = function(x)
{
  is.finite(x) & x > 0
}
isFinitePositiveInteger = function(x)
{
  isFiniteInteger(x) & isFinitePositive(x)
}
isProb = function(x)
{
  is.finite(x) & x >= 0 & x <= 1
}

arrayextend = function(x, along, size, newdimnames)
{
  oldDims = dim(x)
  newdim = oldDims
  newdim[along] = newdim[along] + size
  rnames = dimnames(x)
  rnames[along] = newdimnames
  newarray = tensorA::to.tensor(NA, dim=rnames)
  args = list(newarray)
  args = append(args, lapply(oldDims, function(x) 1:x))
  args[[length(args)+1]] = x
  do.call('[<-',args)
}

#' Slices a tensorA::tensor
#'
#' Slices a tensorA::tensor, preserving the dimnames.
#' This is a workaround for a buggy implementation of
#' [.tensor as of tensorA v0.36.
#' @param x tensor array to be sliced
#' @param ... arguments used to subset the sensor array
sliceTensor <- function(x,...) {
  class(x) <- 'array'
  x = tensorA::as.tensor.default(x[...])
  class(x) = append(class(x),'array')
  x
}


