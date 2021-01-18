#' Simulate false signal probability alpha given control limit for CUSUM charts
#'
#' @export
#' @import checkmate
#' @import stats
#' @param failure_probability Double. Baseline failure probability
#' @param n_patients Integer. Number of patients in monitoring period /sample size
#' @param odds_multiplier Double. Odds multiplier of adverse event under the alternative hypothesis (<1 looks for decreases)
#' @param n_simulation Integer. Number of simulation runs
#' @param limit Double. Control limit for signalling performance change
#' @param seed Integer. Seed for RNG
#' @return Returns False signal probability of specified CUSUM chart.
#' @examples
#'
#' #
#' # control limit can be obtained with cusum_limit_sim(),
#' # here it is set to an arbitrary value (2.96)
#'
#' # simulate false positive probability of CUSUM
#' cusum_alpha_sim(
#'   failure_probability = 0.05,
#'   n_patients = 100,
#'   odds_multiplier = 2,
#'   n_simulation = 10000,
#'   limit = 2.96,
#'   seed = 2046
#' )
cusum_alpha_sim <- function(failure_probability, n_patients, odds_multiplier, n_simulation, limit, seed = NULL) {

  ## Check user input ####
  assert_numeric(failure_probability, lower = 0, upper = 1, finite = TRUE, any.missing = FALSE, len = 1)
  if (failure_probability > 0.5) {
    failure_probability <- 1 - failure_probability
    warning("Accepted failure probability failure_probability will be recoded to 1-failure_probability when > 0.5.")
  }

  assert_integer(as.integer(n_patients), lower = 1, any.missing = FALSE, len = 1)

  assert_numeric(odds_multiplier, lower = 0, finite = TRUE, any.missing = FALSE, len = 1)

  assert_integer(as.integer(n_simulation), lower = 1, any.missing = FALSE, len = 1)
  
  assert_numeric(limit, finite = TRUE, any.missing = FALSE, len = 1)
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
  
  assert_integer(as.integer(seed), lower = 0, upper = Inf, any.missing = TRUE, max.len = 1)
  
  ## Simulate CUSUM Charts ####
  cs_sim <- function(i, npat = n_patients, p = failure_probability, or = odds_multiplier) {
    p.0 <- p
    o.0 <- p.0 / (1 - p.0)
    o.1 <- o.0 * or
    p.1 <- o.1 / (1 + o.1)

    y <- rbinom(npat, 1, p.0)
    w.t <- y * log(p.1 / p.0) + (1 - y) * log((1 - p.1) / (1 - p.0))
    c.t <- vector(mode = "numeric", length = npat)
    if (or > 1){
      c.t[1] <- max(c(0, c.t[1] + w.t[1]))
      for (i in 2:npat) c.t[i] <- max(c(0, c.t[i - 1] + w.t[i]))
      return(max(c.t))
    } else {
      c.t[1] <- min(c(0, c.t[1] - w.t[1]))
      for (i in 2:npat) c.t[i] <- min(c(0, c.t[i - 1] - w.t[i]))
      return(min(c.t))
    }
  }

  suppressWarnings(RNGversion("3.5.0"))

  set.seed(seed)
  rl <- lapply(1:n_simulation, cs_sim)

  ## Estimate false signal probability ####
  x <- ecdf(unlist(rl))
  
  if (odds_multiplier > 1){
    return(1- x(limit))
  } else {
      return(x(limit))
  }

  
}
