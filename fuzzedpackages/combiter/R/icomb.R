#' Combination Iterator
#' @description Create an iterator for all combinations k integers out of 1 through n.
#' @param n positive integer
#' @param k positive integer no greater than n
#' @return iterator object
#' @export
#' @examples
#' x <- icomb(5, 3)
#' ct <- 0
#' while (hasNext(x))
#' {
#'   ct <- ct + 1
#'   i <- nextElem(x)
#'   cat(sprintf("%3d : %s\n", ct, paste0(i, collapse = " ")))
#' }
icomb <- function(n, k)
{
  stopifnot(length(n) == 1L)
  stopifnot(n > 0L)
  stopifnot((n %% 1) == 0)

  stopifnot(length(k) == 1L)
  stopifnot(k > 0L)
  stopifnot((k %% 1) == 0)
  stopifnot(k <= n)
  n <- as.integer(n)
  k <- as.integer(k)


  obj <- recursiveiter(nextFunc = function(i,n,k) NextComb(i, n),
                       prevFunc = function(i,n,k) PrevComb(i, n),
                       firstFunc = function(n,k) 1:k,
                       lastFunc = function(n,k) (n-k+1):n,
                       n = n, k = k)
  class(obj) <- c("icomb", class(obj))
  obj
}



#' @export
#' @rdname icomb
#' @param values iterable (subsettable by \code{[})
#' @details \itemize{
#' \item{\code{icomb} iterates through integer vectors}
#' \item{\code{icombv} iterates through general values}
#' }
#' @examples
#'
#' as.list(icombv(c("A", "G", "C"), 2))
icombv <- function(values, k)
{
  n <- length(values)
  obj <- icomb(n,k)
  nextElem <- function() values[obj$nextElem()]
  prevElem <- function() values[obj$prevElem()]
  hasNext  <- function() obj$hasNext()
  hasPrev  <- function() obj$hasPrev()
  getFirst <- function() values[obj$getFirst()]
  getLast  <- function() values[obj$getLast()]

  out <- list(nextElem=nextElem, prevElem=prevElem,
              hasNext=hasNext, hasPrev=hasPrev,
              getFirst=getFirst, getLast=getLast)
  class(out) <- c("icombv", "recursiveiter", "abstractiter", "iter")
  out
}