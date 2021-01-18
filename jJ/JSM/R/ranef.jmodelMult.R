# Random Effects Estimate for Joint Model with NMRE

ranef.jmodelMult <- function(object, ...) {
  if (!inherits(object, "jmodelMult"))
    stop("Can only be used for 'jmodelMult' objects. \n")
  object$est.bi
}