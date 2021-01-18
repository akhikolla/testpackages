#' @useDynLib eddington
#' @importFrom Rcpp sourceCpp
NULL

# Function definitions ------------------------------------------------------

#' Determine if a dataset satisfies a specified Eddington number
#'
#' Indicates whether a certain Eddington number is satisfied, given the data.
#'
#' @inheritParams E_num
#' @param candidate The Eddington number to test for.
#' @seealso \code{\link{E_cum}}, \code{\link{E_next}}, \code{\link{E_num}},
#'   \code{\link{E_req}}
#' @return A logical vector of length 1.
#' @export
E_sat <- function(rides, candidate) {

  length(rides) >= candidate && sum(rides >= candidate) >= candidate

}

#' Determine the number of additional rides required to achieve a specified
#' Eddington number
#'
#' Determine the number of additional rides required to achieve a specified
#' Eddington number.
#'
#' @inheritParams E_num
#' @param candidate The Eddington number to test for.
#' @seealso \code{\link{E_cum}}, \code{\link{E_next}}, \code{\link{E_num}},
#'   \code{\link{E_sat}}
#' @return An integer vector of length 1. Returns \code{0L} if \emph{E} is
#'   already achieved.
#' @export
E_req <- function(rides, candidate) {

  max(as.integer(candidate) - sum(rides >= candidate), 0L)

}

# Custom print methods -------------------------------------------------------
#' @export
print.E_next <- function(x, ...) {

  out_string <- sprintf(
    "Your current Eddington Number is %i. You need %i %s of %i or greater to get to an Eddington number of %i.",
    x$E,
    x$req,
    if (x$req > 1) "rides" else "ride",
    x$E + 1,
    x$E + 1
  )

  cat(strwrap(out_string), sep = "\n")

}
