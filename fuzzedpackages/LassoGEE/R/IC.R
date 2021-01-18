#' @title Information Criterion for selecting the tuning parameter.
#' @description  Information Criterion for a fitted LassoGEE object with
#'  the AIC, BIC, or GCV criteria.
#' @param obj A fitted LassoGEE object.
#' @param criterion The criterion by which to select the regularization parameter.
#' One of "AIC", "BIC", "GCV", "AICc", or "EBIC"; default is "BIC".
#' @return
#'   \item{IC}{The calculated model selection criteria}
#' @references Gao, X., and Yi, G. Y. (2013). Simultaneous model selection and estimation
#' for mean and association structures with clustered binary data. Stat, 2(1), 102-118.
#' @export
#'
IC <- function(obj, criterion=c("BIC","AIC","GCV","AICc","EBIC")) {
  criterion <- match.arg(criterion)
  scor <- obj$S
  scorv <- obj$Smat
  scorv <- scorv%*%t(scorv)
  quasill <- t(scor)%*%ginv(scorv)%*%scor ## t(scor)%*%scor
  betaest <- obj$betaest
  df <- as.double(sum(betaest!=0))
  p <- length(obj$betaest)
  j <- if(obj$family$family=="gaussian") df - 2 else df - 1
  N <- length(unique(obj$id))
  # N <- obj$nobs

  IC <- switch(criterion,
               AIC = quasill + 2*df,
               BIC = quasill + log(N)*df,
               GCV = (1/N) * (-2) * as.double(quasill) / (1-df/N)^2,
               AICc = quasill + 2*df + 2*df*(df+1)/(N-df-1),
               EBIC = quasill + log(N)*df + 2*(lgamma(p+1) - lgamma(j+1) - lgamma(p-j+1)))

  return(IC=IC)
}

