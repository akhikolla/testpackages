#' Bernstein Polynomial Based Regression Object Summary BPPH Bayes
#'
#' @export
#' @param ... further arguments passed to or from other methods
#' @method print summary.bpph.bayes
#' @return none

print.summary.bpph.bayes <-
  function(...){
  cat("Bayesian Bernstein Polynomial based Proportional Hazards model\n")
  print.summary.spbp.bayes(...)
}
