#' Bernstein Polynomial Based Regression Object Summary BPPO Bayes
#'
#' @export
#' @param ... further arguments passed to or from other methods
#' @method print summary.bppo.bayes
#' @return none

print.summary.bppo.bayes <-
  function(...){
  cat("Bayesian Bernstein Polynomial based Proportional Odds model\n")
  print.summary.spbp.bayes(...)
}
