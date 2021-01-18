
AIC.jmodelTM <-  function (object, ...) {
  - 2 * object$logLik + 2 * nrow(object$Vcov)
}
  
