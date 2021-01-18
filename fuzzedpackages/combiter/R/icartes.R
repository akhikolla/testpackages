#' Cartesian Product Iterator
#'
#' Create an iterator going through Cartesian product of several items.
#' @param nvec integer vector of number of items
#' @return iterator object
#' @export
#' @examples
#' x <- icartes(c(3, 2, 4))
#' ct <- 0
#' while (hasNext(x))
#' {
#'   ct <- ct + 1
#'   i <- nextElem(x)
#'   cat(sprintf("%3d : %s\n", ct, paste0(i, collapse = " ")))
#' }
icartes <- function(nvec)
{
  stopifnot(is.numeric(nvec) > 0)
  stopifnot(length(nvec) > 0)
  stopifnot(all(nvec >= 1))
  stopifnot(all(nvec %% 1 == 0))
  nvec <- as.integer(nvec)

  obj <- recursiveiter(nextFunc = function(i, nvec) NextCartes(i, nvec),
                       prevFunc = function(i, nvec) PrevCartes(i, nvec),
                       firstFunc = function(nvec) rep(1L, length(nvec)),
                       lastFunc = function(nvec) nvec,
                       nvec = nvec)
  class(obj) <- c("icartes", class(obj))
  obj
}



#' @export
#' @rdname icartes
#' @param ... set of iterables (subsettable by \code{[})
#' @details \itemize{
#' \item{\code{icartes} iterates through all combinations of integers}
#' \item{\code{icartesv} iterates through all combinations of general values}
#' }
#' @examples
#'
#' x <- icartesv(Month=c("Jan", "Feb", "Mar"),
#'               Loc=c("NY", "LA"),
#'               By=c("car", "plane", "bus"))
#' as.list(x)
icartesv <- function(...)
{
  value_set <- list(...)
  nvec <- unlist(lapply(value_set, length))
  obj <- icartes(nvec)
  nextElem <- function() Map(`[`, value_set, obj$nextElem())
  prevElem <- function() Map(`[`, value_set, obj$prevElem())
  hasNext  <- function() obj$hasNext()
  hasPrev  <- function() obj$hasPrev()
  getFirst <- function() Map(`[`, value_set, obj$getFirst())
  getLast  <- function() Map(`[`, value_set, obj$getLast())

  out <- list(nextElem=nextElem, prevElem=prevElem,
              hasNext=hasNext, hasPrev=hasPrev,
              getFirst=getFirst, getLast=getLast)
  class(out) <- c("isubsetv", "recursiveiter", "abstractiter", "iter")
  out
}