#' Calculate exact control limit given false signal probability alpha for CUSUM charts for very small sample sizes
#'
#' This function only works for very small sample sizes (<= 15),
#' as it permutes through all possible outcome sequences and estimates the percentage of runs that reach a specific CUSUM values.
#'
#' @export
#' @import checkmate
#' @param failure_probability Double. Baseline failure probability
#' @param n_patients Integer. Number of patients in monitoring period /sample size
#' @param odds_multiplier Double. Odds multiplier of adverse event under the alternative hypothesis (<1 looks for decreases)
#' @param alpha Double. False signal probability of CUSUM
#' @return Returns the control limit for signalling performance change for small sample sizes (double)
#' @examples
#'
#' # calculate exact control limits for alpha = 0.05
#' cusum_limit_exact(
#'   failure_probability = 0.1,
#'   n_patients = 10,
#'   odds_multiplier = 2,
#'   alpha = 0.05
#' )
cusum_limit_exact <- function(n_patients,
                              failure_probability,
                              odds_multiplier,
                              alpha) {
  #
  # This function calculates the exact distribution of the CUSUM
  # for a hospital with n_patients patients, the in control failure probability failure_probability
  # and the smallest inacceptable failure probability pA
  #

  assert_numeric(failure_probability, lower = 0, upper = 1, finite = TRUE, any.missing = FALSE, len = 1)
  if (failure_probability > 0.5) {
    failure_probability <- 1 - failure_probability
    warning("Baseline failure probability failure_probability will be recoded to 1 - failure_probability when > 0.5.")
  }

  assert_integer(as.integer(n_patients), lower = 1, any.missing = FALSE, len = 1)
  if (n_patients > 20) {
    stop("Exact calculation only works for very small sample sizes of around <= 10 (check ?cusum_limit_exact for more information). \nPlease abort if calculation takes to long. ")
  } else if (n_patients > 12) {
    message("Exact calculation only works for very small sample sizes of around <= 10 (check ?cusum_limit_exact for more information). \nPlease abort if calculation takes to long. ")
  }

  assert_numeric(odds_multiplier, lower = 1, finite = TRUE, any.missing = FALSE, len = 1)

  if (odds_multiplier == 1) {
    stop("CUSUM detects no process change (odds_multiplier = 1).")
  }

  assert_numeric(alpha, lower = 0, upper = 1, finite = TRUE, any.missing = FALSE, len = 1)

  p.0 <- failure_probability
  o.0 <- p.0 / (1 - p.0)
  o.1 <- o.0 * odds_multiplier
  p.1 <- o.1 / (1 + o.1)

  outcome <- make_all_outcomes(npat_outcome = n_patients)

  p_failure <- apply(outcome, 1, function(kk, pp) {
    return(prod(ifelse(kk == 1, pp, 1 - pp)))
  }, pp = failure_probability)

  which_rfc <- function(xx, mm) {
    res <- apply(mm, 2, function(yy, cc) {
      return(sum(yy >= cc) > 0)
    }, cc = xx)
    
    return(res)
  }
  cs <- apply(outcome, 1, calc_cusum, c0 = failure_probability, cA = p.1)
  limit <- unique(as.vector(cs))
  which.res <- lapply(limit, which_rfc, mm = cs)

  cs_distr <- lapply(which.res, function(yy, pp) {
    return(sum(pp[yy]))
  }, pp = p_failure)

  res <- cbind(limit, unlist(cs_distr))
  res <- res[sort.list(res[, 1]), ]
  return(res[which.min(abs(alpha - res[, 2])), 1])
}

#' Make all outcomes
#'
#' creates all possible sequences of outcomes for a sample size
#'
#' @param npat_outcome Number of patients (sample sizes)
#' @return Returns matrix of possible sequences
make_all_outcomes <- function(npat_outcome) {

  #
  # 	This function creates all possible sequences of outcomes
  #

  m <- matrix(0:1, ncol = 1)

  for (ii in 2:npat_outcome) {
    m <- cbind(rbind(m, m), rep(0:1, c(1, 1) * 2^(ii - 1)))
  }

  return(m)
}

#' Calculate CUSUM
#'
#' This function calculates the CUSUM chart for the given sequence of successes and failures
#'
#' @param x vector of outcomes
#' @param c0 accepted failure probability
#' @param cA smallest detectable failure probability
#' @return Returns matrix of possible sequences
calc_cusum <- function(x, c0, cA) {
  wt <- ifelse(x == 0, log((1 - cA) / (1 - c0)), log(cA / c0))

  j <- length(wt)
  ct <- rep(NA, j)
  ct[1] <- max(c(0, wt[1]))

  for (ii in 2:j) {
    ct[ii] <- max(c(0, ct[ii - 1] + wt[ii]))
  }
  return(ct)
}


