
#' Implementation of Math group generics for lvec
#'
#' @param x an object of type \code{\link{lvec}}.
#' @param ... passed on to the corresponding R functions
#'
#' @details
#' Math is group generic implementing the following functions: 
#' \code{\link{abs}}, \code{\link{sign}}, \code{\link{sqrt}},
#' \code{\link{floor}}, \code{\link{ceiling}}, \code{\link{trunc}},
#' \code{\link{round}}, \code{\link{signif}} \code{\link{exp}}, 
#' \code{\link{log}}, \code{\link{expm1}}, \code{\link{log1p}},
#' \code{\link{cos}}, \code{\link{sin}}, \code{\link{tan}},
#' \code{\link{cospi}}, \code{\link{sinpi}}, \code{\link{tanpi}},
#' \code{\link{acos}}, \code{\link{asin}}, \code{\link{atan}},
#' \code{\link{cosh}}, \code{\link{sinh}}, \code{\link{tanh}},
#' \code{\link{acosh}}, \code{\link{asinh}}, \code{\link{atanh}},
#' \code{\link{lgamma}}, \code{\link{gamma}}, \code{\link{digamma}}, 
#' \code{\link{trigamma}}, \code{\link{cumsum}}, \code{\link{cumprod}}, 
#' \code{\link{cummax}}, \code{\link{cummin}}. For more information
#' see \code{\link{Math}}.
#'
#' @return
#' Returns an \code{link{lvec}} of the same length as the input.
#'
#' @export
Math.lvec <- function(x, ...) {
  fun <- function(d, ...) {
    do.call(.Generic, list(x = d, ...))
  }
  cumulative_functions <- c("cumsum", "cumprod", "cummax", "cummin")
  if (.Generic %in% cumulative_functions) {
    elementwise_cumultative(x, fun, ...)
  } else {
    elementwise(x, fun, ...)
  }
}



#' @import lvec
elementwise_cumultative <- function(x, cumfun, ...) {
  chunks <- chunk(x)
  result <- NULL
  last_element <- NULL

  for (c in chunks) {
    d <- as_rvec(lget(x, range = c))
    r <- cumfun(c(last_element, d), ...)
    if (!is.null(last_element)) r <- r[-1]
    last_element <- utils::tail(r, 1)
    r <- as_lvec(r)
    if (is.null(result)) {
      result <- r
      length(result) <- length(x)
    } else {
      if (lvec_type(r) == 'character' && strlen(r) > strlen(result)) {
        warning("Changing maximum string length")
        strlen(result) <- strlen(r)
      }
      lset(result, range = c, r)
    }
  }
  result
}

