
BIC.jmodelTM <-  function (object, ...) {
  - 2 * object$logLik + log(object$n) * nrow(object$Vcov)
}