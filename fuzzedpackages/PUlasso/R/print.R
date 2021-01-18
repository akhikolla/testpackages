#' @export
#' @method print PUfit
print.PUfit<-function(x,...){
  cat("\nCall: ", deparse(x$call), "\n")
  cat("\nOptimization Method: ", deparse(x$optResult$method), "\n")
  cat("\nMax nIters: ", deparse(x$optResult$maxit), "\n\n")
}

#' @export
#' @method print cvPUfit
print.cvPUfit<-function(x,...){
  cat("\nCall: ", deparse(x$call), "\n")
  cat("\nOptimization Method: ", deparse(x$PUfit$optResult$method), "\n")
  cat("\nMax nIters: ", deparse(x$PUfit$optResult$maxit), "\n\n")
}
