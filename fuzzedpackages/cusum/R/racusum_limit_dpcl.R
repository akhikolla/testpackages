#' Dynamic Probability Control Limits (DPCL)
#'
#'
#' Set DPCL for risk-adjusted Bernoulli CUSUM Charts
#' @export
#' @import stats
#' @import data.table
#' @param patient_risks Double. Vector of patient risk scores (individual risk of adverse event)
#' @param N Integer. Number of simulation runs
#' @param odds_multiplier Double. Odds multiplier of adverse event under the alternative hypothesis (<1 looks for decreases)
#' @param alpha Double. False signal probability of RA-CUSUM
#' @param seed Integer. Seed for RNG
#' @references Zhang, Xiang & Woodall, William. (2016). Dynamic Probability Control Limits for Lower and Two-Sided Risk-Adjusted Bernoulli CUSUM Charts. Quality and Reliability Engineering International. 10.1002/qre.2044.
#' @return Returns vector of dynamic control limit for signalling performance change (double)
#'
#' @examples
#' patient_risks <- runif(100, min = 0.1, max = 0.8)
#'
#' dpcl <- racusum_limit_dpcl(
#'   patient_risks = patient_risks,
#'   N = 1000,
#'   odds_multiplier = 2,
#'   alpha = 0.05,
#'   seed = 32423
#' )
#' 
#' plot(dpcl, type = "l")
racusum_limit_dpcl <- function(patient_risks, N = 100000, odds_multiplier = 2, alpha, seed = NULL) {
  ## Check user input ####
  assert_numeric(patient_risks, lower = 0, upper = 1, finite = TRUE, any.missing = FALSE, min.len = 1)

  assert_integer(as.integer(N), lower = 1, any.missing = FALSE, len = 1)
  
  assert_numeric(odds_multiplier, lower = 0, finite = TRUE, any.missing = FALSE, len = 1)
  if (odds_multiplier < 1) {
    message("CUSUM is set to detect process improvements (odds_multiplier < 1). ")
  }
  if (odds_multiplier == 1) {
    warning("CUSUM is set to detect no process change (odds_multiplier = 1).")
  }

  assert_numeric(alpha, lower = 0, upper = 1, finite = TRUE, any.missing = FALSE, len = 1)
  
  assert_integer(as.integer(seed), lower = 0, upper = Inf, any.missing = TRUE, max.len = 1)
  
  cs <- 0

  M <- (N * (1 - alpha))

  h <- vector(length = length(patient_risks))

  set.seed(seed)
  for (i in 1:length(patient_risks)) {
    pi <- patient_risks[i]
    yi <- rbinom(N, 1, pi)
    ws <- log((1 - pi + 1 * pi) / (1 - pi + odds_multiplier * pi)) # survival
    wf <- log(((1 - pi + 1 * pi) * odds_multiplier) / ((1 - pi + odds_multiplier * pi) * 1)) # failure
    w <- ifelse(yi == 1, wf, ws)

    cs_i <- NULL
    for (j in 1:N) {
      ct_i <- sample(cs, 1)
      cs_i[j] <- max(0, ct_i + w[j])
    }
    cs <- sort(cs_i)
    h[i] <- cs[M]
    cs <- cs[cs <= h[i]]
  }


  return(h)
}
