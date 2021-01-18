
#' Factory of Iterators defined by Recursive Transition Functions
#'
#' @description This is a constructor for custom iterator objects.
#' It requires four functions, "next", "prev", "first", and "last", and
#' additional parameters.
#'
#' The state of the constructor is characterized by the variable \code{i}.
#' The "next" and "prev" function must take \code{i} and the parameters
#' and return the next and previous state variables respectively.  The behavior where there is no more state left is arbitrary.
#'
#' The "first" and "last" functions must take the additional parameters and
#' return the initial and last state variables respectively.
#'
#' The created object is an iterator of class \code{recursiveiter}, which inherits
#' \code{abstractiter} and \code{iter}.
#' It can be used with \code{\link[foreach]{foreach}} and accepts \code{\link{as.list}} conversion.
#'
#' @param nextFunc,prevFunc Functions that take the iterator state and the parameters \code{...} and returns the next or previous state
#' @param firstFunc,lastFunc  Functions that take the parameters \code{...}
#' and returns the first or last state of the iteration
#' @param ... additional parameters of the iterator
#' @export
#' @examples
#' fibiter <- recursiveiter(
#'   nextFunc = function(i) if (length(i)==1 && i==0) 1 else
#'                          if (length(i)==1 && i==1) c(1,1) else
#'                          c(sum(i), i[1]),
#'   prevFunc = NULL, firstFunc = function() 0, lastFunc = function() Inf)
#' for (k in 1:20) cat(nextElem(fibiter)[1], "")
#' @return iterator object
recursiveiter <- function(nextFunc, prevFunc, firstFunc, lastFunc, ...)
{
  i <- NULL

  .hasNext <- function() { is.null(i) || !identical(i, lastFunc(...)) }
  .hasPrev <- function() { is.null(i) || !identical(i, firstFunc(...)) }

  .nextElem <- function()
  {
    if (is.null(i)) {
      i <<- firstFunc(...)
    } else if (!.hasNext()) {
      stop("StopIteration")
    } else {
      i <<- nextFunc(i, ...)
    }
    return(i)
  }

  .prevElem <- function()
  {
    if (is.null(i)) {
      i <<- lastFunc(...)
    } else if (!.hasPrev()) {
      stop("StopIteration")
    } else {
      i <<- prevFunc(i, ...)
    }
    return(i)
  }

  .getFirst <- function() firstFunc(...)
  .getLast  <- function() lastFunc(...)

  #self <- environment()
  #class(self) <- "combinatiter"
  #return(self)
  out <- list(nextElem=.nextElem, prevElem=.prevElem,
              hasNext=.hasNext, hasPrev=.hasPrev,
              getFirst=.getFirst, getLast=.getLast)
  class(out) <- c("recursiveiter", "abstractiter", "iter")
  out
}


## S3 methods

#' @export
nextElem.recursiveiter <- function(obj, ...) { obj$nextElem() }

#' @export
prevElem.recursiveiter <- function(obj, ...) { obj$prevElem() }

#' @export
hasNext.recursiveiter <- function(obj, ...) { obj$hasNext() }

#' @export
hasPrev.recursiveiter <- function(obj, ...) { obj$hasPrev() }

#' @export
getFirst.recursiveiter <- function(obj, ...) { obj$getFirst() }

#' @export
getLast.recursiveiter <- function(obj, ...) { obj$getLast() }
