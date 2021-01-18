#' Permutation Iterator
#' @description Create an iterator for all permutations of size k of integers 1 to n.
#' @param n positive integer
#' @param k positive integer
#' @return iterator object
#' @export
#' @examples
#' x <- iperm(3)
#' ct <- 0
#' while (hasNext(x))
#' {
#'   ct <- ct + 1
#'   i <- nextElem(x)
#'   cat(sprintf("%3d : %s\n", ct, paste0(i, collapse = " ")))
#' }
#'
iperm <- function(n, k = n)
{
  stopifnot(length(n) == 1L)
  stopifnot(n > 0L)
  stopifnot((n %% 1) == 0)
  stopifnot(length(k) == 1L)
  stopifnot(k > 0L)
  stopifnot((k %% 1) == 0)
  stopifnot(n >= k)
  k <- as.integer(k)
  n <- as.integer(n)

  obj <- recursiveiter(nextFunc = function(i, n, k) NextPerm(i, n),
                       prevFunc = function(i, n, k) PrevPerm(i, n),
                       firstFunc = function(n, k) 1:k,
                       lastFunc = function(n, k) n:(n-k+1),
                       n = n, k = k)
  class(obj) <- c("iperm", class(obj))
  obj
}



#' @rdname iperm
#' @export
#' @param values iterable (subsettable by \code{[})
#' @details \itemize{
#' \item{\code{iperm} iterates through integer vectors}
#' \item{\code{ipermv} iterates through general values}
#' }
#' @examples
#'
#' as.list(ipermv(c("R", "G", "B")))
ipermv <- function(values, k = length(values))
{
  n <- length(values)
  obj <- iperm(n,k)
  nextElem <- function() values[obj$nextElem()]
  prevElem <- function() values[obj$prevElem()]
  hasNext  <- function() obj$hasNext()
  hasPrev  <- function() obj$hasPrev()
  getFirst <- function() values[obj$getFirst()]
  getLast  <- function() values[obj$getLast()]

  out <- list(nextElem=nextElem, prevElem=prevElem,
              hasNext=hasNext, hasPrev=hasPrev,
              getFirst=getFirst, getLast=getLast)
  class(out) <- c("ipermv", "recursiveiter", "abstractiter", "iter")
  out
}