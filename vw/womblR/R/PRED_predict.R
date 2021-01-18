###Function to get model fit diagnostics given a STBDwDM object
#'
#' predict.STBDwDM
#'
#' Predicts future observations from the \code{\link{STBDwDM}} model.
#'
#' @param object a \code{\link{STBDwDM}} model object for which predictions
#'  are desired from.
#'
#' @param NewTimes a numeric vector including desired time(s) points for prediction.
#'
#' @param ... other arguments.
#'
#' @details \code{predict.STBDwDM} uses Bayesian krigging to predict vectors at future
#'  time points. The function returns the krigged observed outcomes along with the
#'  observational level parameters (\code{mu}, \code{tau}, and \code{alpha}).
#'
#' @return \code{predict.STBDwDM} returns a list containing the following objects.
#'
#'   \describe{
#'
#'   \item{\code{MuTauAlpha}}{A \code{list} containing three matrices, \code{mu},
#'   \code{tau} and \code{alpha}. Each matrix is dimension \code{NKeep x s}, where
#'   \code{s} is the number of new time points. Each matrix contains posterior
#'   samples obtained by Bayesian krigging.}
#'
#'   \item{\code{Y}}{A \code{list} containing \code{s} posterior predictive distribution
#'   matrices. Each matrix is dimension \code{NKeep x s}, where \code{s}
#'   is the number of new time points. Each matrix is obtained through Bayesian krigging.}
#'
#'   }
#'
#' @author Samuel I. Berchuck
#' @export
###Prediction function for spBDwDM function
predict.STBDwDM <- function(object, NewTimes, ...) {

  ###Check Inputs
  if (missing(object)) stop('"object" is missing')
  if (!is.STBDwDM(object)) stop('"object" must be of class STBDwDM')
  if (missing(NewTimes)) stop('"NewTimes" is missing')
  if (!is.numeric(NewTimes)) stop('NewTimes must be a vector')
  if (any(is.na(NewTimes))) stop("NewTimes may have no missing values")
  if (any(!is.finite(NewTimes))) stop("NewTimes must have strictly finite entries")
  if (!all(NewTimes >= 0)) stop('NewTimes vector has at least one negative entry')

  ###Set seed for reproducibility
  set.seed(54)

  ###Set data objects
  DatObj <- object$datobj
  Nu <- DatObj$Nu
  M <- DatObj$M

  ###Create updated distance matrix
  TimeFixed <- DatObj$Time
  Time <- sort(c(TimeFixed, NewTimes))
  TimeDist <- abs(outer(Time, Time, "-" ))
  NNewVisits <- length(NewTimes)
  NewVisits <- OriginalVisits <- NULL
  for (i in 1:NNewVisits) NewVisits <- c(NewVisits, which(NewTimes[i] == Time) - 1)
  for (i in 1:Nu) OriginalVisits <- c(OriginalVisits, which(TimeFixed[i] == Time) - 1)

  ###Update DatObj
  DatObj$NewVisits <- NewVisits
  DatObj$OriginalVisits <- OriginalVisits
  DatObj$TimeDist <- TimeDist
  DatObj$NNewVisits <- NNewVisits

  ###Set mcmc object
  NKeep <- dim(object$phi)[1]

  ###Create parameter object
  Para <- list()
  Para$Mu <- object$mu
  Para$Tau2 <- object$tau2
  Para$Alpha <- object$alpha
  Para$Delta <- object$delta
  Para$T <- object$T
  Para$Phi <- object$phi

  ###Obtain samples of mu, tau and alpha using Bayesian krigging
  ThetaKrig <- ThetaKrigging(DatObj, Para, NKeep)

  ###Obtain samples of observed Y
  YKrig <- YKrigging(DatObj, ThetaKrig, NKeep)

  ###Format theta samples for output
  FutureThetaArray <- array(ThetaKrig, dim = c(NKeep, 3, NNewVisits))
  Mu <- matrix(FutureThetaArray[ , 1 , ], ncol = NNewVisits)
  Tau2 <- matrix(exp(FutureThetaArray[ , 2 , ])^2, ncol = NNewVisits)
  Alpha <- matrix(exp(FutureThetaArray[ , 3 , ]), ncol = NNewVisits)
  colnames(Mu) <- paste0("mu", NewVisits + 1)
  colnames(Tau2) <- paste0("tau2", NewVisits + 1)
  colnames(Alpha) <- paste0("alpha", NewVisits + 1)
  LevelOneKrig <- list(mu = Mu, tau2 = Tau2, alpha = Alpha)

  ###Format Y samples for output
  FutureYArray <- array(t(YKrig), dim = c(NKeep, M, NNewVisits))
  PpdKrig <- list()
  for (i in 1:NNewVisits) PpdKrig[[i]] <- FutureYArray[ , , i]
  PpdKrig <- lapply(PpdKrig, f <- function(x) {colnames(x) <- paste0("loc", 1:M); return(x);})
  names(PpdKrig) <- paste0("y", NewVisits + 1)

  ###Return formated samples
  return(list(MuTauAlpha = LevelOneKrig, Y = PpdKrig))

}
