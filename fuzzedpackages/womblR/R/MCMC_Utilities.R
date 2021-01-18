###Function for summarizing the raw MCMC samples-------------------------------------------------------------------
FormatSamples <- function(DatObj, RawSamples) {

  ###Set data objects
  Nu <- DatObj$Nu

  ###Format raw samples
  Mu <- t(RawSamples[1 : Nu, ])
  Tau2 <- t(RawSamples[(1 : Nu) + Nu, ])
  Alpha <- t(RawSamples[(1 : Nu) + 2 * Nu, ])
  Delta <- t(RawSamples[(3 * Nu + 1) : (3 * Nu + 3), ])
  T <- t(RawSamples[(3 * Nu + 4) : (3 * Nu + 9), ])
  Phi <- t(RawSamples[(3 * Nu + 10), ,drop = FALSE])
  colnames(Mu) <- paste0("mu", 1 : Nu)
  colnames(Tau2) <- paste0("tau2", 1 : Nu)
  colnames(Alpha) <- paste0("alpha", 1 : Nu)
  colnames(Delta) <- paste0("delta", 1 : 3)
  colnames(T) <- c("t11", "t21", "t22", "t31", "t32", "t33")
  colnames(Phi) <- "phi"
  Out <- list(Mu = Mu, Tau2 = Tau2, Alpha = Alpha, Delta = Delta, T = T, Phi = Phi)
  return(Out)

}



###Function for creating a data object that contains objects needed for ModelFit-----------------------------------
OutputDatObj <- function(DatObj) {

  ###Collect needed objects
  DatObjOut <- list(M = DatObj$M,
                    Nu = DatObj$Nu,
                    Z = DatObj$Z,
                    AdjacentEdgesBoolean = DatObj$AdjacentEdgesBoolean,
                    W = DatObj$W,
                    EyeM = DatObj$EyeM,
                    OneM = DatObj$OneM,
                    OneNu = DatObj$OneNu,
                    YStarWide = DatObj$YStarWide,
                    Rho = DatObj$Rho,
                    FamilyInd = DatObj$FamilyInd,
                    ScaleY = DatObj$ScaleY,
                    YObserved = DatObj$YObserved,
                    ScaleDM = DatObj$ScaleDM,
                    TempCorInd = DatObj$TempCorInd,
                    WeightsInd = DatObj$WeightsInd,
                    Time = DatObj$Time)
  return(DatObjOut)

}



###Function for creating a data augmentation object that contains objects needed for ModelFit----------------------
OutputDatAug <- function(DatAug) {

  ###Collect needed objects
  DatAugOut <- list(NBelow = DatAug$NBelow,
                    NBelowCount = DatAug$NBelowCount,
                    TobitIndeces = DatAug$TobitIndeces,
                    YStarNonZero = DatAug$YStarNonZero)
  return(DatAugOut)

}



###Function for summarizing Metropolis objects post sampler--------------------------------------------------------
SummarizeMetropolis <- function(DatObj, MetrObj, McmcObj) {

  ###Set data object
  Nu <- DatObj$Nu

  ###Set MCMC object
  NSims <- McmcObj$NSims

  ###Set Metropolis objects
  MetropTheta2 <- MetrObj$MetropTheta2
  AcceptanceTheta2 <- MetrObj$AcceptanceTheta2
  MetropTheta3 <- MetrObj$MetropTheta3
  AcceptanceTheta3 <- MetrObj$AcceptanceTheta3
  MetropPhi <- MetrObj$MetropPhi
  AcceptancePhi <- MetrObj$AcceptancePhi

  ###Summarize and output
  TuningParameters <- c(MetropTheta2, MetropTheta3, MetropPhi)
  AcceptanceCount <- c(AcceptanceTheta2, AcceptanceTheta3, AcceptancePhi)
  AcceptancePcts <- AcceptanceCount / NSims
  MetrSummary <- cbind(AcceptancePcts, TuningParameters)
  rownames(MetrSummary) <- c(paste0("Theta2", 1 : Nu), paste0("Theta3", 1 : Nu), "Phi")
  colnames(MetrSummary) <- c("acceptance", "tuner")
  return(MetrSummary)

}


###Verify the class of our regression object------------------------------------------------------------------------
#' is.STBDwDM
#'
#' \code{is.STBDwDM} is a general test of an object being interpretable as a
#' \code{\link{STBDwDM}} object.
#'
#' @param x object to be tested.
#'
#' @details The \code{\link{STBDwDM}} class is defined as the regression object that
#'  results from the \code{\link{STBDwDM}} regression function.
#'
#' @export
is.STBDwDM <- function(x) {
  identical(attributes(x)$class, "STBDwDM")
}


