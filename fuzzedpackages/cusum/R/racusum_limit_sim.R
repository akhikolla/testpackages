#' Simulate control limit given false signal probability alpha for RA-CUSUM charts
#'
#' False-signal-probability-simulation of Control Limits h for risk-adjusted CUSUM charts
#' @export
#' @import checkmate
#' @import stats
#' @param patient_risks Double. Vector of patient risk scores (individual risk of adverse event)
#' @param n_simulation Integer. Number of simulation runs
#' @param odds_multiplier Double. Odds multiplier of adverse event under the alternative hypothesis (<1 looks for decreases)
#' @param alpha Double. False signal probability of RA-CUSUM
#' @param seed Integer. Seed for RNG
#' @return Returns the control limit for signalling performance change (double)
#' @examples
#'
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
#' # simulate control limits for alpha = 0.05
#' racusum_limit_sim(patient_risks,
#'   odds_multiplier = 2,
#'   n_simulation = 1000,
#'   alpha = 0.05,
#'   seed = 2046
#' )
racusum_limit_sim <- function(patient_risks, odds_multiplier, n_simulation, alpha, seed = NULL) {

  ## Check user input ####
  assert_numeric(patient_risks, lower = 0, upper = 1, finite = TRUE, any.missing = FALSE, min.len = 1)

  assert_numeric(odds_multiplier, lower = 0, finite = TRUE, any.missing = FALSE, len = 1)
  if (odds_multiplier < 1) {
   # message("CUSUM detects process improvements (odds_multiplier < 1). ")
  }
  if (odds_multiplier == 1) {
    stop("CUSUM detects no process change (odds_multiplier = 1).")
  }

  assert_integer(as.integer(n_simulation), lower = 1, any.missing = FALSE, len = 1)
  
  assert_numeric(alpha, lower = 0, upper = 1, finite = TRUE, any.missing = FALSE, len = 1)

  assert_integer(as.integer(seed), lower = 0, upper = Inf, any.missing = TRUE, max.len = 1)
  
  ## Calculate risk distribution ####
  n <- length(patient_risks)

  freq_score <- table(patient_risks)
  freq_score <- as.data.frame(freq_score)


  freq <- round(as.numeric(freq_score$Freq / sum(freq_score$Freq)), digits = 1000)
  risk <- sort(unique(patient_risks)) # corresponding risk scores


  if (sum(freq) != 1) {
    warning("Sum of risk frequence != 1.")
  }

  ## Simulate CUSUM runs ####
  cs_sim <- function(i) {
    p.0 <- sample(risk, size = n, replace = T, prob = freq)


    y <- rbinom(n, 1, p.0)
    ws <- log(1 / (1 + (odds_multiplier - 1) * p.0)) # success (non death) case
    wf <- log(odds_multiplier / (1 + (odds_multiplier - 1) * p.0)) # failure (death) case

    w.t <- ifelse(y == 1, wf, ws)
    c.t <- vector(mode = "numeric", length = n)
    
    if (odds_multiplier > 1){
      c.t[1] <- max(c(0, c.t[1] + w.t[1]))
      for (i in 2:n) c.t[i] <- max(c(0, c.t[i - 1] + w.t[i]))
      return(max(c.t))
    } else {
      c.t[1] <- min(c(0, c.t[1] - w.t[1]))
      for (i in 2:n) c.t[i] <- min(c(0, c.t[i - 1] - w.t[i]))
      return(min(c.t))
    }
  }
  suppressWarnings(RNGversion("3.5.0"))

  set.seed(seed)
  rl <- lapply(1:n_simulation, cs_sim)

  ## Estimate Alpha ####
  if (odds_multiplier > 1){
    q <- quantile(unlist(rl), 1 - alpha)
    
  } else {
    q <- quantile(unlist(rl), alpha)
    
  }

  return(as.numeric(q))
}
