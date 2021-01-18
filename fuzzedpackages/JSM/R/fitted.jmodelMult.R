# Fitted Values for Joint Model with NMRE

fitted.jmodelMult <- function (object, process = c("Longitudinal", "Survival"),
                               type = c("Marginal", "Conditional"), ...) {
  if (!inherits(object, "jmodelMult"))
    stop("Can only be used for 'jmodelMult' objects. \n")
  process <- match.arg(process)
  if (process == "Survival")
    stop("fitted() is not implemented for the survival process yet. \n")
  type <- match.arg(type)
  fitY <- as.vector(object$dataMat$B %*% object$coefficients$gamma)
  names(fitY) <- object$dataMat$IDName
  if (type == "Marginal") fitY else if (type == "Conditional") fitY * object$est.bi[object$dataMat$ID, ]
}