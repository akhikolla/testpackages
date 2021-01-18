## generic functions
### nextElem ... from iterators
### hasNext  ... from itertools
### prevElem, hasPrev, getFirst, getLast ... defined here

#' @importFrom iterators nextElem
#' @export
iterators::nextElem

#' @importFrom itertools hasNext
#' @export
itertools::hasNext

#' Get Previous Element of Iterator
#' @description \code{prevElem} is a generic funcion to move an
#' iterator object one step backward.
#' @param obj an R object
#' @param ... additional arguments
#' @return iterator value
#' @export
prevElem <- function(obj, ...)
{
  UseMethod("prevElem")
}


#' Does This Iterator Have A Previous Element
#' @description \code{hasPrev} is a generic function that indicates if the
#' iterator has another element backward.
#' @param obj an R object
#' @param ... additional arguments
#' @return Logical value indicating whether the iterator has a previous element.
#' @export
hasPrev <- function(obj, ...)
{
  UseMethod("hasPrev")
}


#' First Value of Iterator
#' @description \code{getFirst} is a generic function that returns the
#' first value of iterators
#' @param obj an R object
#' @param ... additional arguments
#' @return iterator value, format dependes on the objects
#' @export
getFirst <- function(obj, ...)
{
  UseMethod("getFirst")
}

#' Last Value of Iterator
#' @description \code{getFirst} is a generic function that returns the
#' last value of iterators
#' @param obj an R object
#' @param ... additional arguments
#' @return iterator value, format dependes on the objects
#' @export
getLast <- function(obj, ...)
{
  UseMethod("getLast")
}

