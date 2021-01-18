# Fitted Values for Joint Model with LME

fitted.jmodelTM <- function (object, process = c("Longitudinal", "Survival"), 
                             type = c("Marginal", "Conditional"), ...) {
  if (!inherits(object, "jmodelTM"))
    stop("Can only be used for 'jmodelTM' objects. \n")
  process <- match.arg(process)
  if (process == "Survival")
    stop("fitted() is not implemented for the survival process yet. \n")
  type <- match.arg(type)
  fitY <- as.vector(object$dataMat$X %*% object$coefficients$beta)
  names(fitY) <- object$dataMat$IDName
  Zb <- as.vector(rowSums(object$dataMat$Z * object$est.bi[object$dataMat$ID, ]))
  if (type == "Marginal") fitY else if (type == "Conditional") fitY + Zb
}