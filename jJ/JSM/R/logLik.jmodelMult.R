
logLik.jmodelMult <-  function (object, ...) {
  if (!inherits(object, "jmodelMult"))
    stop("Only used for 'jmodelMult' objects.\n")
  out <- object$logLik
  attr(out, "df") <- nrow(object$Vcov)
  attr(out, "n") <- object$n
  class(out) <- "logLik"
  out
}
