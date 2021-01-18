#' Subset Iterator
#' @description Create an iterator for all subsets of integers 1 through n.
#' @param n positive integer
#' @return iterator object
#' @export
#' @examples
#' x <- isubset(3)
#' ct <- 0
#' while (hasNext(x))
#' {
#'   ct <- ct + 1
#'   i <- nextElem(x)
#'   cat(sprintf("%3d : %s\n", ct, paste0(i, collapse = " ")))
#' }
isubset <- function(n)
{
  stopifnot(length(n) == 1L)
  stopifnot(n > 0L)
  stopifnot((n %% 1) == 0)
  n <- as.integer(n)

  obj <- recursiveiter(nextFunc = function(i, n) NextSubset(i, n),
                       prevFunc = function(i, n) PrevSubset(i, n),
                       firstFunc = function(n) integer(0),
                       lastFunc = function(n) 1:n,
                       n = n)
  class(obj) <- c("isubset", class(obj))
  obj
}


#' @export
#' @rdname isubset
#' @param values iterable (subsettable by \code{[})
#' @details \itemize{
#' \item{\code{isubset} iterates through integer vectors}
#' \item{\code{isubsetv} iterates through general values}
#' }
#' @examples
#'
#' as.list(isubsetv(letters[1:4]))
isubsetv <- function(values)
{
  n <- length(values)
  obj <- isubset(n)
  nextElem <- function() values[obj$nextElem()]
  prevElem <- function() values[obj$prevElem()]
  hasNext  <- function() obj$hasNext()
  hasPrev  <- function() obj$hasPrev()
  getFirst <- function() values[obj$getFirst()]
  getLast  <- function() values[obj$getLast()]

  out <- list(nextElem=nextElem, prevElem=prevElem,
              hasNext=hasNext, hasPrev=hasPrev,
              getFirst=getFirst, getLast=getLast)
  class(out) <- c("isubsetv", "recursiveiter", "abstractiter", "iter")
  out
}

