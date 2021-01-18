#' Weights for observations
#'
#' Calculate standard CUSUM weights
#' @export 
#' @param patient_outcomes Integer. Vector of binary patient outcomes (0,1)
#' @param probability_ae Double. Baseline failure probability for adverse event in non-risk-adjusted case, vector of patient risk scores for risk-adjustment.
#' @param odds_multiplier Double. Odds multiplier of adverse event under the alternative hypothesis (<1 looks for decreases)
#' 

weights_t <- function(patient_outcomes,
                      probability_ae,
                      odds_multiplier = 2) {
  assert_numeric(odds_multiplier, lower = 0, finite = TRUE, any.missing = FALSE, len = 1)
  assert_integer(patient_outcomes, lower = 0, upper = 1, any.missing = FALSE, min.len = 1)
  
  if (length(probability_ae) == 1) {
    assert_numeric(probability_ae, lower = 0, upper = 1, finite = TRUE, any.missing = FALSE, len = 1)
    if (probability_ae > 0.5) {
      probability_ae <- 1 - probability_ae
    }
    ## non-risk-adjusted

    o0 <- probability_ae / (1 - probability_ae)
    oA <- odds_multiplier * o0
    pA <- oA / (1 + oA)

    wf <- log(pA / probability_ae)
    ws <- log((1 - pA) / (1 - probability_ae))
  } else {
    assert_numeric(probability_ae, lower = 0, upper = 1, finite = TRUE, any.missing = FALSE, min.len = 1)
    ## risk-adjusted
    ws <- log((1)/ (1 - probability_ae + odds_multiplier * probability_ae))
    wf <- log((odds_multiplier) / ((1 - probability_ae + odds_multiplier * probability_ae) * 1))
  }
  
  w <- ifelse(patient_outcomes == 1, wf, ws)
  return(w)
}
