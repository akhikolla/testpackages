#' Non-risk-adjusted CUSUM Charts
#'
#' Calculate non-risk-adjusted CUSUM charts for performance data
#' @export
#' @import checkmate
#' @import stats
#' @import graphics
#' @param failure_probability Double. Baseline failure probability
#' @param patient_outcomes Integer. Vector of binary patient outcomes (0,1)
#' @param limit Double. Control limit for signalling performance change
#' @param weights Double. Optional vector of weights, if empty, standard CUSUM weights are calculated with weights_t
#' @param odds_multiplier Double. Odds multiplier of adverse event under the alternative hypothesis (<1 looks for decreases)
#' @param reset Logical. Reset the CUSUM after a signal to 0; defaults to TRUE
#' @examples
#'
#' # control limit can be obtained with cusum_limit_sim(),
#' # here it is set to an arbitrary value (2.96)
#'
#' # CUSUM of in-control process
#' # simulate patient outcomes
#' set.seed(2046)
#' patient_outcomes <- as.logical(rbinom(n = 100, size = 1, prob = 0.05))
#'
#'
#' cs_ic <- cusum(
#'   failure_probability = 0.05,
#'   patient_outcomes,
#'   limit = 2.96
#' )
#'
#' # CUSUM of out-of-control process
#' # simulate patient outcome
#' set.seed(2046)
#' patient_outcomes <- as.logical(rbinom(n = 100, size = 1, prob = 0.2))
#'
#' cs_oc <- cusum(
#'   failure_probability = 0.05,
#'   patient_outcomes,
#'   limit = 2.96
#' )
cusum <- function(failure_probability, patient_outcomes, limit, weights = NULL, odds_multiplier = 2, reset = TRUE) {

  ## Check user input ####
  assert_numeric(failure_probability, lower = 0, upper = 1, finite = TRUE, any.missing = FALSE, len = 1)
  if (failure_probability > 0.5) {
    failure_probability <- 1 - failure_probability
    warning("Accepted failure probability failure_probability will be recoded to 1-failure_probability when > 0.5.")
  }

  patient_outcomes <- as.integer(patient_outcomes)
  assert_integer(patient_outcomes, lower = 0, upper = 1, any.missing = FALSE, min.len = 1)

  assert_numeric(limit, finite = TRUE, any.missing = FALSE, len = 1)

  if (length(weights) > 0){
    assert_numeric(weights, lower = -1, upper = 1, finite = TRUE, any.missing = FALSE, min.len = 1)
    if (length(weights) != length(patient_outcomes)) {
      stop("Length weights and patient outcomes of unequal size.")
    }
  }
  
  assert_numeric(odds_multiplier, lower = 0, finite = TRUE, any.missing = FALSE, len = 1)
  if (odds_multiplier < 1) {
   # message("CUSUM is set to detect process improvements (odds_multiplier < 1). ")

    if (limit > 0) {
      warning("Control limit should be negative to signal process improvements.")
    }
  }
  if (odds_multiplier == 1) {
    stop("CUSUM is set to detect no process change (odds_multiplier = 1).")
  } 
  if (odds_multiplier > 1){
    if (limit < 0) {
      warning("Control limit should be positive to signal process deteriorations.")
    }
  }

  assert_logical(reset, any.missing = FALSE, len = 1)

  ## Calculate CUSUM Chart ####
  npat <- length(patient_outcomes)

  ct <- 0
  cs <- matrix(0, nrow = npat, ncol = 5)
  

  if (length(weights) == 0){
    w <-  weights_t(patient_outcomes,
                    probability_ae = failure_probability,
                    odds_multiplier)
  } else {
    w <- weights
  }


  for (ii in 1:npat) {
    
    if (odds_multiplier > 1) {
      ct <- max(0, ct + w[ii])
      if (ct >= limit) {
        cs[ii, 4] <- 1 # store signal
        if (reset == TRUE) ct <- 0
      } else {
        cs[ii, 4] <- 0
      }
    } else if (odds_multiplier < 1) {
      ct <- min(0, ct - w[ii])
      if (ct <= limit) {
        # test for signal
        cs[ii, 4] <- 1 # store signal

        if (reset == TRUE) ct <- 0
      } else {
        cs[ii, 4] <- 0
      }
    }
    
    cs[ii, 1] <- ii # store patient id
    cs[ii, 2] <- failure_probability
    cs[ii, 3] <- ct # store CUSUM value
  }

  cs[, 5] <- limit

  cs <- as.data.frame(cs)
  names(cs) <- c("t", "failure_probability", "ct", "signal", "limit")

  class(cs) <- c("cusum", "data.frame")

  return(cs)
}
