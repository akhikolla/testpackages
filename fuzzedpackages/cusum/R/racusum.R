#' Risk-adjusted CUSUM Charts
#'
#' Calculate risk-adjusted CUSUM charts for performance data
#' @export
#' @import checkmate
#' @import stats
#' @import graphics
#' @param patient_risks Double. Vector of patient risk scores (individual risk of adverse event)
#' @param patient_outcomes Integer. Vector of binary patient outcomes (0,1)
#' @param limit Double. Control limit for signalling performance change
#' @param weights Double. Optional vector of weights, if empty, standard CUSUM weights are calculated with weights_t
#' @param odds_multiplier Double. Odds multiplier of adverse event under the alternative hypothesis (<1 looks for decreases)
#' @param reset Logical. Reset the CUSUM after a signal to 0; defaults to TRUE
#' @param limit_method "constant" or "dynamic"
#' @examples
#' # Patients risks are usually known from Phase I.
#' # If not, these risk scores can be simulated.
#'
#' # define possible patient risk scores
#' risks <- c(0.001, 0.01, 0.1, 0.002, 0.02, 0.2)
#'
#' # sample risk population of size n = 100
#' set.seed(2046)
#' patient_risks <- sample(x = risks, size = 100, replace = TRUE)
#'
#' # control limit can be obtained with racusum_limit_sim(),
#' # here it is set to an arbitrary value (2.96),
#' # or dynamic control limits with racusum_limit_dpcl()
#'
#' ##### RA-CUSUM of in-control process
#' # simulate patient outcome for performace as expected
#' set.seed(2046)
#' patient_outcomes <- as.logical(rbinom(
#'   n = 100,
#'   size = 1,
#'   prob = patient_risks
#' ))
#'
#' racusum(patient_risks,
#'   patient_outcomes,
#'   limit = 2.96
#' )
#'
#' #### RA-CUSUM of out-of-control process
#' # simulate patient outcome for deviating performance
#'
#' set.seed(2046)
#' patient_outcomes <- as.logical(rbinom(n = 100, size = 1, prob = patient_risks * 2))
#' #'
#' racusum(patient_risks,
#'   patient_outcomes,
#'   limit = 2.96
#' )
racusum <- function(patient_risks, patient_outcomes, limit, weights = NULL, odds_multiplier = 2, reset = TRUE, limit_method = c("constant", "dynamic")) {
  npat <- length(patient_risks)
  
  ## Check user input ####
  assert_numeric(patient_risks, lower = 0, upper = 1, finite = TRUE, any.missing = FALSE, min.len = 1)
  if (length(patient_risks) != length(patient_outcomes)) {
    stop("Length patient_risks and patient_outcomes of unequal size.")
  }
  
  if (length(weights) > 0){
    assert_numeric(weights, lower = -1, upper = 1, finite = TRUE, any.missing = FALSE, min.len = 1)
    if (length(weights) != length(patient_outcomes)) {
      stop("Length weights and patient outcomes of unequal size.")
    }
  }
  
  patient_outcomes <- as.integer(patient_outcomes)
  assert_integer(patient_outcomes, lower = 0, upper = 1, any.missing = FALSE, min.len = 1)
  
  if (!missing(limit_method)) {
    warning("argument limit_method is deprecated and not needed anymore.",
            call. = FALSE
    )
  }
  
  if (length(limit) == 1) {
    limit <- rep(limit, length.out = npat)
  }
  
  
  assert_numeric(odds_multiplier, lower = 0, finite = TRUE, any.missing = FALSE, len = 1)
  if (odds_multiplier < 1) {
    # message("CUSUM is set to detect process improvements (odds_multiplier < 1). ")
    
    if (mean(limit) > 0) {
      warning("Control limit should be negative to signal process improvements.")
    }
  }
  if (odds_multiplier == 1) {
    # warning("CUSUM is set to detect no process change (odds_multiplier = 1).")
  }
  if (odds_multiplier > 1){
    if (mean(limit) < 0) {
      warning("Control limit should be positive to signal process deteriorations.")
    }
  }
  
  
  assert_logical(reset, any.missing = FALSE, len = 1)
  
  ## Calculate RA-CUSUM Chart ####
  npat <- length(patient_risks)
  
  if (length(weights) == 0){
    w <- weights_t(patient_outcomes,
                   probability_ae = patient_risks,
                   odds_multiplier)
  } else {
    w <- weights
  }
  
  ct <- 0 # initial CUSUM value
  
  cs <- matrix(0, nrow = npat, ncol = 5) # storage matrix for alarms
  # p <- patient_risks
  # 
  # 
  # # CUSUM weights
  # ws <- log((1 - p + 1 * p) / (1 - p + odds_multiplier * p))
  # wf <- log(((1 - p + 1 * p) * odds_multiplier) / ((1 - p + odds_multiplier * p) * 1))
  # 
  # w <- ifelse(patient_outcomes == 1, wf, ws) # weights based on outcome

  for (ii in 1:npat) {
    if (odds_multiplier > 1) {
      ct <- max(0, ct + w[ii])

      if (ct >= limit[ii]) {
        # test for signal
        cs[ii, 4] <- 1 # store signal
        cs[, 5] <- limit[ii]
        if (reset == 1) ct <- 0
      } else {
        cs[ii, 4] <- 0
      }
    } else if (odds_multiplier < 1) {
      ct <- min(0, ct - w[ii])

      if (ct <= limit[ii]) {
        # test for signal
        cs[ii, 4] <- 1 # store signal
        cs[, 5] <- limit[ii]
        if (reset == 1) ct <- 0
      } else {
        cs[ii, 4] <- 0
      }
    }
    cs[ii, 1] <- ii # store patient id
    cs[ii, 2] <- patient_risks[ii] # store patient risk
    cs[ii, 3] <- ct # store CUSUM value
  }



  cs <- as.data.frame(cs)
  names(cs) <- c("t", "p", "ct", "signal", "limit")

  class(cs) <- c("cusum", "data.frame")


  return(cs)
}
