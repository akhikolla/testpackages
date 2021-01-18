###Function for summarizing the raw MCMC samples-------------------------------------------------------------------
FormatSamples <- function(DatObj, RawSamples) {

  ###Set data objects
  M <- DatObj$M
  PhiIndeces <- DatObj$PhiIndeces
  tNu <- DatObj$tNu
  t1 <- DatObj$t1

  ###Format raw samples
  RawSamples <- t(RawSamples)
  Alpha <- RawSamples[, 1, drop = FALSE]
  Delta <- RawSamples[, 2:6]
  Sigma <- RawSamples[, 7:21]
  Phi <- RawSamples[, 22:((5 * M) + 21)]
  Beta0 <- Phi[, PhiIndeces[1, ] + 1]
  Beta1 <- Phi[, PhiIndeces[2, ] + 1]
  Lambda0 <- Phi[, PhiIndeces[3, ] + 1]
  Lambda1 <- Phi[, PhiIndeces[4, ] + 1]
  Eta <- Phi[, PhiIndeces[5, ] + 1]
  Theta <- apply(Eta, 2, function(x) pmax(pmin(x, tNu), t1))
  colnames(Alpha) <- "Alpha"
  colnames(Delta) <- paste0("Delta", 1:5)
  colnames(Sigma) <- c(paste0("Sigma", 1:5, "1"), paste0("Sigma", 2:5, "2"), paste0("Sigma", 3:5, "3"), paste0("Sigma", 4:5, "4"), "Sigma55")
  colnames(Beta0) <- paste0("Beta0(", 1:M, ")")
  colnames(Beta1) <- paste0("Beta1(", 1:M, ")")
  colnames(Lambda0) <- paste0("Lambda0(", 1:M, ")")
  colnames(Lambda1) <- paste0("Lambda1(", 1:M, ")")
  colnames(Eta) <- paste0("Eta(", 1:M, ")")
  colnames(Theta) <- paste0("Theta(", 1:M, ")")
  Out <- list(Alpha = Alpha, Delta = Delta, Sigma = Sigma, Beta0 = Beta0, Beta1 = Beta1, Lambda0 = Lambda0, Lambda1 = Lambda1, Eta = Eta, Theta = Theta)
  return(Out)

}



###Function for creating a data object that contains objects needed for ModelFit-----------------------------------
OutputDatObj <- function(DatObj) {

  ###Collect needed objects
  DatObjOut <- list(M = DatObj$M,
                    Nu = DatObj$Nu,
                    AdjacentEdgesBoolean = DatObj$AdjacentEdgesBoolean,
                    W = DatObj$W,
                    EyeM = DatObj$EyeM,
                    EyeNu = DatObj$EyeNu,
                    OneM = DatObj$OneM,
                    OneN = DatObj$OneN,
                    OneNu = DatObj$OneNu,
                    YStarWide = DatObj$YStarWide,
                    Rho = DatObj$Rho,
                    FamilyInd = DatObj$FamilyInd,
                    ScaleY = DatObj$ScaleY,
                    YObserved = DatObj$YObserved,
                    ScaleDM = DatObj$ScaleDM,
                    Time = DatObj$Time,
                    TimeVec = DatObj$TimeVec,
                    YObserved = DatObj$YObserved,
                    tNu = DatObj$tNu,
                    t1 = DatObj$t1,
                    XThetaInd = DatObj$XThetaInd,
                    N = DatObj$N,
                    EyeN = DatObj$EyeN)
  return(DatObjOut)

}



###Function for creating a data augmentation object that contains objects needed for ModelFit----------------------
OutputDatAug <- function(DatAug) {

  ###Collect needed objects
  DatAugOut <- list(NBelow = DatAug$NBelow,
                    NBelowList = DatAug$NBelowList,
                    TobitBooleanMat = DatAug$TobitBooleanMat,
                    YStarNonZeroList = DatAug$YStarNonZero)
  return(DatAugOut)

}



###Function for summarizing Metropolis objects post sampler--------------------------------------------------------
SummarizeMetropolis <- function(DatObj, MetrObj, MetropRcpp, McmcObj) {

  ###Set data object
  M <- DatObj$M

  ###Set MCMC object
  NSims <- McmcObj$NSims

  ###Set Metropolis objects
  MetropLambda0Vec <- MetropRcpp$MetropLambda0Vec
  AcceptanceLambda0Vec <- MetropRcpp$AcceptanceLambda0Vec
  MetropLambda1Vec <- MetropRcpp$MetropLambda1Vec
  AcceptanceLambda1Vec <- MetropRcpp$AcceptanceLambda1Vec
  MetropEtaVec <- MetropRcpp$MetropEtaVec
  AcceptanceEtaVec <- MetropRcpp$AcceptanceEtaVec
  MetropAlpha <- MetropRcpp$MetropAlpha
  AcceptanceAlpha <- MetropRcpp$AcceptanceAlpha
  OriginalTuners <- MetrObj$OriginalTuners

  ###Summarize and output
  TuningParameters <- c(MetropLambda0Vec, MetropLambda1Vec, MetropEtaVec, MetropAlpha)
  AcceptanceCount <- c(AcceptanceLambda0Vec, AcceptanceLambda1Vec, AcceptanceEtaVec, AcceptanceAlpha)
  AcceptancePcts <- AcceptanceCount / NSims
  MetrSummary <- cbind(AcceptancePcts, TuningParameters, OriginalTuners)
  rownames(MetrSummary) <- c(paste0("Lambda0(", 1:M, ")"), paste0("Lambda1(", 1:M, ")"), paste0("Eta(", 1:M, ")"), "Alpha")
  colnames(MetrSummary) <- c("Acceptance", "PilotAdaptedTuners", "OriginalTuners")
  return(MetrSummary)

}



###Verify the class of our regression object------------------------------------------------------------------------
#' is.spCP
#'
#' \code{is.spCP} is a general test of an object being interpretable as a
#' \code{\link{spCP}} object.
#'
#' @param x object to be tested.
#'
#' @details The \code{\link{spCP}} class is defined as the regression object that
#'  results from the \code{\link{spCP}} regression function.
#'
#' @export
is.spCP <- function(x) {
  identical(attributes(x)$class, "spCP")
}