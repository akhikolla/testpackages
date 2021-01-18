# Residuals for Joint Model with NMRE

residuals.jmodelMult <- function(object, process = c("Longitudinal", "Survival"),
                                 type = c("Marginal", "Conditional", "Standardized-Marginal", "Standardized-Conditional"), ...) {
  if (!inherits(object, "jmodelMult"))
    stop("Can only be used for 'jmodelMult' objects. \n")
  process <- match.arg(process)
  if (process == "Survival")
    stop("residuals() is not implemented for the survival process yet. \n")
  type <- match.arg(type)
  Y <- object$dataMat$Y
  B <- object$dataMat$B
  ID <- object$dataMat$ID
  ncb <- ncol(B)
  
  fitY <- as.vector(B %*% object$coefficients$gamma)
  fittedVal <- if (type %in% c("Marginal", "Standardized-Marginal")) {
    fitY
  } else {
    fitY * object$est.bi[ID, ]
  }
  
  Bsigma <- object$coefficients$Bsigma
  Ysigma <- object$coefficients$Ysigma
  if (type %in% c("Marginal", "Conditional")) {
    residVal <- as.vector(Y - fittedVal)
  } else if (type == "Standardized-Conditional") {
    residVal <- as.vector(Y - fittedVal) / Ysigma
  } else { # type == "Standardized-Marginal"
    residVal <- as.vector((Y - fittedVal) / sqrt(fitY ^ 2 * Bsigma ^ 2 + Ysigma ^2))
  }
  names(residVal) <- object$dataMat$IDName
  residVal
}
