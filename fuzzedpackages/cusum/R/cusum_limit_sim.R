#' Simulate control limit given false signal probability alpha for CUSUM charts
#'
#' @export
#' @import checkmate
#' @import stats
#' @param failure_probability Double. Baseline failure probability
#' @param n_patients Integer. Number of patients in monitoring period /sample size
#' @param odds_multiplier Double. Odds multiplier of adverse event under the alternative hypothesis (<1 looks for decreases)
#' @param n_simulation Integer. Number of simulation runs
#' @param alpha Double. False signal probability of CUSUM
#' @param seed Integer. Seed for RNG
#' @return Returns the control limit for signalling performance change (double)
#' @examples
#'
#' # simulate control limits for alpha = 0.05
#' cusum_limit_sim(
#'   failure_probability = 0.05,
#'   n_patients = 100,
#'   odds_multiplier = 2,
#'   n_simulation = 1000,
#'   alpha = 0.05,
#'   seed = 2046
#' )
cusum_limit_sim <- function(failure_probability, n_patients, odds_multiplier, n_simulation, alpha, seed = NULL) {
  
  ## Check user input ####
  assert_numeric(failure_probability, lower = 0, upper = 1, finite = TRUE, any.missing = FALSE, len = 1)
  if (failure_probability > 0.5) {
    failure_probability <- 1 - failure_probability
    warning("Accepted failure probability failure_probability will be recoded to 1-failure_probability when > 0.5.")
  }
  
  n <- as.integer(n_patients)
  assert_integer(as.integer(n_patients), lower = 1, any.missing = FALSE, len = 1)

  assert_numeric(odds_multiplier, lower = 0, finite = TRUE, any.missing = FALSE, len = 1)
  if (odds_multiplier < 1) {
    #message("CUSUM detects process improvements (odds_multiplier < 1). ")
  }
  if (odds_multiplier == 1) {
    stop("CUSUM detects no process change (odds_multiplier = 1).")
  }

  assert_integer(as.integer(n_simulation), lower = 1, any.missing = FALSE, len = 1)
  
  assert_numeric(alpha, lower = 0, upper = 1, finite = TRUE, any.missing = FALSE, len = 1)

  assert_integer(as.integer(seed), lower = 0, upper = Inf, any.missing = TRUE, max.len = 1)

  ## Simulate CUSUM runs ####
  cs_sim <- function(i, npat = n, p = failure_probability, or = odds_multiplier) {
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

  ## Estimate Alpha ####
  if (odds_multiplier > 1){
    cl <- quantile(unlist(rl), 1 - alpha)
    
  } else {
    cl <- quantile(unlist(rl), alpha)
  }

  return(as.numeric(cl))
}
