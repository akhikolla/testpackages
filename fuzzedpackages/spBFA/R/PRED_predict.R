###Function to get model fit diagnostics given a spBFA object
#'
#' predict.spBFA
#'
#' Predicts future observations from the \code{\link{spBFA}} model.
#'
#' @param object A \code{\link{spBFA}} model object for which predictions
#'  are desired from.
#'
#' @param NewTimes A numeric vector including desired time(s) points for prediction.
#' 
#' @param NewX A matrix including covariates at times \code{NewTimes} for prediction. 
#'  \code{NewX} must have dimension \code{(M x O x NNewVistis) x P}. Where \code{NNewVisits} is the number of temporal 
#'  locations being predicted. The default sets \code{NewX} to \code{NULL}, which assumes that the covariates for all predictions 
#'  are the same as the final time point.
#'  
#' @param NewTrials An array indicating the trials for categorical predictions. The array must have dimension \code{M x C x NNewVisits}
#'  and contain only non-negative integers. The default sets \code{NewTrials} to \code{NULL}, which assumes the trials for all predictions
#'  are the same as the final time point.
#'  
#' @param type A character string indicating the type of prediction, choices include "temporal" and "spatial". Spatial prediction has not been implemented yet.
#'
#' @param Verbose A boolean logical indicating whether progress should be output.
#'
#' @param ... other arguments.
#'
#' @details \code{predict.spBFA} uses Bayesian krigging to predict vectors at future
#'  time points. The function returns the krigged factors (\code{Eta}) and also the observed outcomes (\code{Y}).
#'
#' @return \code{predict.spBFA} returns a list containing the following objects.
#'
#'   \describe{
#'
#'   \item{\code{Eta}}{A \code{list} containing \code{NNewVistis} matrices, one for each new time prediction. Each matrix is dimension \code{NKeep x K}, where
#'   \code{K} is the number of latent factors Each matrix contains posterior samples obtained by Bayesian krigging.}
#'
#'   \item{\code{Y}}{A \code{list} containing \code{NNewVistis} posterior predictive distribution
#'   matrices. Each matrix is dimension \code{NKeep x (M * O)}, where \code{M} is the number of spatial locations and \code{O} the number of observation types.
#'   Each matrix is obtained through Bayesian krigging.}
#'
#'   }
#'
#' @examples
#' ###Load pre-computed regression results
#' data(reg.bfa_sp)
#' 
#' ###Compute predictions
#' pred <- predict(reg.bfa_sp, NewTimes = 3)
#' pred.observations <- pred$Y$Y10 # observed data predictions
#' pred.krig <- pred$Eta$Eta10 # krigged parameters
#' 
#' @author Samuel I. Berchuck
#' @export
###Prediction function for spBFA function
predict.spBFA <- function(object, NewTimes, NewX = NULL, NewTrials = NULL, Verbose = TRUE, type = "temporal", ...) {

  ###Check Inputs
  if (missing(object)) stop('"object" is missing')
  if (!is.spBFA(object)) stop('"object" must be of class spBFA')
  if (missing(NewTimes)) stop('"NewTimes" is missing')
  if (!is.numeric(NewTimes)) stop('NewTimes must be a vector')
  if (any(is.na(NewTimes))) stop("NewTimes may have no missing values")
  if (any(!is.finite(NewTimes))) stop("NewTimes must have strictly finite entries")
  if (!all(NewTimes >= 0)) stop('NewTimes vector has at least one negative entry')
  if (!is.character(type)) stop('"type" must be a character string')
  if (!(type %in% c("temporal", "spatial"))) stop('"type" must be one of "spatial" or "temporal"')
  if (!is.logical(Verbose)) stop('"Verbose" must be a logical')
  
  ###Set seed for reproducibility
  set.seed(54)

  ###Set data objects
  DatObj <- object$datobj
  Nu <- DatObj$Nu
  M <- DatObj$M
  O <- DatObj$O
  P <- DatObj$P
  K <- DatObj$K

  ###Create updated distance matrix
  TimeFixed <- DatObj$Time
  Time <- sort(c(TimeFixed, NewTimes))
  TimeDist <- abs(outer(Time, Time, "-" ))
  NNewVisits <- length(NewTimes)
  NewVisits <- OriginalVisits <- NULL
  for (i in 1:NNewVisits) NewVisits <- c(NewVisits, which(NewTimes[i] == Time) - 1)
  for (i in 1:Nu) OriginalVisits <- c(OriginalVisits, which(TimeFixed[i] == Time) - 1)

  ###Get covariates
  if (is.null(NewX)) {
    XNu <- object$datobj$X[object$datobj$Indeces == (object$datobj$Nu - 1), ]
    NewX <- do.call("rbind", rep(list(XNu), NNewVisits))
  } else {
    if (!is.matrix(NewX)) stop('NewX must be a matrix')
    if (dim(NewX)[1] != (M * O * NNewVisits)) stop("NewX: Must be a matrix with dimension (M x O x NNewVisits) x P")
    if (dim(NewX)[2] != P) stop("NewX: Must be a matrix with dimension (M x O x NNewVisits) x P")
    if (any(is.na(NewX))) stop("NewX may have no missing values")
    if (any(!is.finite(NewX))) stop("NewX must have strictly finite entries")
    NewX <- NewX
  }
  
  ###Update DatObj
  DatObj$NewVisits <- NewVisits
  DatObj$OriginalVisits <- OriginalVisits
  DatObj$TimeDist <- TimeDist
  DatObj$NNewVisits <- NNewVisits
  DatObj$EyeK <- diag(DatObj$K)
  DatObj$NewX <- NewX
  
  ###Create Trials object
  if (DatObj$C == 0) {
    DatObj$Trials <- array(0, dim = c(M, O, NNewVisits))
  }
  if ((DatObj$C > 0) & is.null(NewTrials)) {
    Trials <- array(dim = c(M, DatObj$C, NNewVisits))
    for (n in 1:NNewVisits) Trials[, , n] <- DatObj$Trials[, , Nu]
  }
  if ((DatObj$C > 0) & !is.null(NewTrials)) {
    Trials <- NewTrials
    if (!is.array(Trials)) stop('Trials must be an array')
    if (dim(Trials)[1] != M) stop("NewTrials: Must be an array with dimension M x C x NNewVisits")
    if (dim(Trials)[2] != DatObj$C) stop("NewTrials: Must be an array with dimension M x C x NNewVisits")
    if (dim(Trials)[3] != NNewVisits) stop("NewTrials: Must be an array with dimension M x C x NNewVisits")
    if (any(is.na(Trials))) stop("Trials may have no missing values")
    if (any(!is.finite(Trials))) stop("Trials must have strictly finite entries")
    if (!isTRUE(all(Trials == floor(Trials)))) stop("Trials must have integers only")
    if (any(Trials < 0)) stop("Trials must contain non-negative integers only")
  }
  
  ###Set mcmc object
  NKeep <- dim(object$rho)[1]

  ###Create parameter object
  Para <- list()
  Para$Psi <- object$psi
  Para$Upsilon <- object$upsilon
  Para$Lambda <- object$lambda
  Para$Eta <- object$eta
  if (DatObj$P > 0) Para$Beta <- object$beta
  if (DatObj$P == 0) Para$Beta <- matrix(0, ncol = 0, nrow = NKeep)
  if (is.null(object$sigma2)) Para$Sigma2 <- matrix(1)
  if (!is.null(object$sigma2)) Para$Sigma2 <- object$sigma2
  
  ###Obtain samples of eta using Bayesian krigging
  EtaKrig <- EtaKrigging(DatObj, Para, NKeep, Verbose)

  ###Obtain samples of observed Y
  YKrig <- YKrigging(DatObj, Para, EtaKrig, NKeep, Verbose)

  ###Format theta samples for output
  EtaOut <- list()
  Eta <- array(t(EtaKrig), dim = c(K, NKeep, NNewVisits))
  for (n in 1:NNewVisits) EtaOut[[n]] <- t(Eta[, , n])
  for (n in 1:NNewVisits) colnames(EtaOut[[n]]) <- paste0("Eta", 1:K, "_", Nu + n)
  names(EtaOut) <- paste0("Eta", Nu + 1:NNewVisits)
  
  ###Format Y samples for output
  YOut <- list()
  YInd <- expand.grid(1:M, 1:O)
  for (n in 1:NNewVisits) YOut[[n]] <- t(YKrig[, n, ])
  for (n in 1:NNewVisits) colnames(YOut[[n]]) <- paste0("Y_", Nu + n, "_", YInd[, 2], "_", YInd[, 1])
  names(YOut) <- paste0("Y", Nu + 1:NNewVisits)  
    
  ###Return formated samples
  return(list(Eta = EtaOut, Y = YOut))

}
