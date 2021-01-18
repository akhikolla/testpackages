# Residuals for Joint Model with LME

residuals.jmodelTM <- function(object, process = c("Longitudinal", "Survival"),
                               type = c("Marginal", "Conditional", "Standardized-Marginal", "Standardized-Conditional"), ...) {
  if (!inherits(object, "jmodelTM"))
    stop("Can only be used for 'jmodelTM' objects. \n")
  process <- match.arg(process)
  if (process == "Survival")
    stop("residuals() is not implemented for the survival process yet. \n")
  type <- match.arg(type)
  Y <- object$dataMat$Y
  X <- object$dataMat$X
  Z <- object$dataMat$Z
  ID <- object$dataMat$ID
  ncx <- ncol(X)
  ncz <- ncol(Z)
  
  fitY <- as.vector(X %*% object$coefficients$beta)
  Zb <- as.vector(rowSums(Z * object$est.bi[ID, ]))
  fittedVal <- if (type %in% c("Marginal", "Standardized-Marginal")) {
    fitY
  } else {
    fitY + Zb
  }

  if (type %in% c("Marginal", "Conditional")) {
    residVal <- as.vector(Y - fittedVal)
  } else if (type == "Standardized-Conditional") {
    residVal <- as.vector(Y - fittedVal) / object$coefficients$Ysigma
  } else { # type == "Standardized-Marginal"
    Bsigma <- object$coefficients$Bsigma
    Ysigma <- object$coefficients$Ysigma
    residVal <- unlist(lapply(split(cbind(Z, as.vector(Y - fittedVal)), ID), function (x) {
      tempMat <- matrix(x, ncol = (ncz + 1))
      Zi <- tempMat[, 1:ncz]
      ri <- tempMat[, (ncz + 1)]
      Vi <- Zi %*% Bsigma %*% t(Zi) 
      diag(Vi) <- diag(Vi) + Ysigma ^ 2
      solve(chol(Vi)) %*% ri
    }))
  }
  names(residVal) <- object$dataMat$IDName
  residVal
}
