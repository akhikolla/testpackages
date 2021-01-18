#' Bernstein Polynomial Based Regression Object Summary BPPH MLE
#'
#' @export
#' @param ... further arguments passed to or from other methods
#' @method print summary.bpph.mle
#' @return none

print.summary.bpph.mle <-
  function(...){
  cat("Bernstein Polynomial based Proportional Hazards model\n")
  print.summary.spbp.mle(...)
}
