#####################################################
# Name: credibleIntervals.R                         #
# Author: Poppy Miller <p.miller@lancaster.ac.uk> & #
# Chris Jewell <c.jewell@lancaster.ac.uk>           #
# Created: 20161206                                 #
# Copyright: Poppy Miller & Chris Jewell 2016       #
# Purpose: Implements credible interval calculation #
#####################################################
#' @importFrom SPIn SPIn

check_alphavalue_CI <- function(alpha) {
  assertthat::assert_that(is.atomic(alpha),
                        is.numeric(alpha),
                        alpha >=0,
                        alpha <= 1)
}

ci_chenShao = function(x, alpha) {
  if (! check_alphavalue_CI(alpha)) stop("alpha is not a single positive numeric value between 0 and 1.")
  n <- length(x)
  sorted <- sort(x)

  upper_pos <- round(n * (1 - alpha))

  ci.lower <- 0
  ci.upper <- 0
  region_size <- Inf

  ## checks all intervals that are n*alpha apart, and chooses the shortest one.
  for (i in 1:(n - upper_pos)) {
    test_interval <- sorted[upper_pos + i] - sorted[i]
    if (test_interval < region_size) {
      region_size <- test_interval
      ci.lower <- sorted[i]
      ci.upper <- sorted[upper_pos + i]
    }
  }

  return(c(
    lower = ci.lower,
    median = stats::median(sorted),
    upper = ci.upper
  ))
}


ci_percentiles <- function(x, alpha) {
  if (! check_alphavalue_CI(alpha)) stop("alpha is not a single positive numeric value between 0 and 1.")
  n <- length(x)
  sorted <- sort(x)

  upper_pos <- round(n * (1 - (alpha / 2)))
  lower_pos <- round(n * (alpha / 2))

  return(c(
    lower = sorted[lower_pos],
    median = stats::median(x),
    upper = sorted[upper_pos]
  ))
}


ci_SPIn <- function(x, alpha) {
  if (! check_alphavalue_CI(alpha)) stop("alpha is not a single positive numeric value between 0 and 1.")
  region <- tryCatch({
    SPIn::SPIn(x, conf = 1 - alpha)$spin
  },
  error = function(cond) {
    print("Error calculating SPIn interval.")
    return(c(NA, NA))
  })
  return(c(
    lower = region[1],
    median = stats::median(x),
    upper = region[2]
  ))
}
