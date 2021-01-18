# Random Effects Estimate for Joint Model with LME

ranef.jmodelTM <- function(object, ...) {
  if (!inherits(object, "jmodelTM"))
    stop("Can only be used for 'jmodelTM' objects. \n")
  object$est.bi
}