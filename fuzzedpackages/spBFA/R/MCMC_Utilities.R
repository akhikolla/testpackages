###Function for summarizing the raw MCMC samples-------------------------------------------------------------------
FormatSamples <- function(DatObj, RawSamples) {

  ###Set data objects
  M <- DatObj$M
  K <- DatObj$K
  Nu <- DatObj$Nu
  O <- DatObj$O
  C <- DatObj$C
  P <- DatObj$P
  GS <- DatObj$GS
  CL <- DatObj$CL
  
  ###Format raw samples
  RawSamples <- t(RawSamples)
  Lambda <- RawSamples[, 1:(O * M * K)]
  Eta <- RawSamples[, (O * M * K + 1):(O * M * K + K * Nu), drop = FALSE]
  if (C == O) Sigma2 <- NULL
  if (C != O) Sigma2 <- RawSamples[, (O * M * K + K * Nu + 1):(O * M * K + K * Nu + M * (O - C))]
  Kappa <- RawSamples[, (O * M * K + K * Nu + M * (O - C) + 1):(O * M * K + K * Nu + M * (O - C) + (O * (O + 1)) / 2), drop = FALSE]
  Delta <- RawSamples[, (O * M * K + K * Nu + M * (O - C) + (O * (O + 1)) / 2 + 1):(O * M * K + K * Nu + M * (O - C) + (O * (O + 1)) / 2 + K), drop = FALSE]
  if (GS == 1) Tau <- matrix(t(apply(Delta, 1, cumprod)), ncol = K)
  if (GS == 0) Tau <- Delta
  Upsilon <- RawSamples[, (O * M * K + K * Nu + M * (O - C) + (O * (O + 1)) / 2 + K + 1):(O * M * K + K * Nu + M * (O - C) + (O * (O + 1)) / 2 + K + (K * (K + 1)) / 2), drop = FALSE]
  Psi <- RawSamples[, (O * M * K + K * Nu + M * (O - C) + (O * (O + 1)) / 2 + K + (K * (K + 1)) / 2 + 1), drop = FALSE]
  Xi <- RawSamples[, (O * M * K + K * Nu + M  * (O - C)+ (O * (O + 1)) / 2 + K + (K * (K + 1)) / 2 + 1 + 1):(O * M * K + K * Nu + M * (O - C) + (O * (O + 1)) / 2 + K + (K * (K + 1)) / 2 + 1 + O * M * K)]
  Rho <- RawSamples[, (O * M * K + K * Nu + M * (O - C) + (O * (O + 1)) / 2 + K + (K * (K + 1)) / 2 + 1 + O * M * K + 1), drop = FALSE]
  if (P > 1) Beta <- RawSamples[, (O * M * K + K * Nu + M * (O - C) + (O * (O + 1)) / 2 + K + (K * (K + 1)) / 2 + 1 + O * M * K + 1 + 1):(O * M * K + K * Nu + M * (O - C) + (O * (O + 1)) / 2 + K + (K * (K + 1)) / 2 + 1 + O * M * K + 1 + P), drop = FALSE]
  LambdaInd <- expand.grid(1:K, 1:M, 1:O)
  colnames(Lambda) <- paste0("Lambda_", LambdaInd[, 3], "_", LambdaInd[, 2], "_", LambdaInd[, 1])
  EtaInd <- expand.grid(1:K, 1:Nu)
  colnames(Eta) <- paste0("Eta", EtaInd[, 2], "_", EtaInd[, 1])
  if (C != O) Sigma2Ind <- expand.grid(which(DatObj$FamilyInd != 3), 1:M)
  if (C != O) colnames(Sigma2) <- paste0("Sigma2_", Sigma2Ind[, 1], "_", Sigma2Ind[, 2])
  KappaInd <- which(lower.tri(apply(matrix(1:O, ncol = 1), 1, function(x) paste0(paste0("Kappa", 1:O, "_"), x)), diag = TRUE), arr.ind = TRUE)
  if (O == 1) colnames(Kappa) <- apply(matrix(1:O, ncol = 1), 1, function(x) paste0(paste0("Kappa", 1:O, "_"), x))[KappaInd[order(KappaInd[, 1]), ]][1]
  if (O > 1) colnames(Kappa) <- apply(matrix(1:O, ncol = 1), 1, function(x) paste0(paste0("Kappa", 1:O, "_"), x))[KappaInd[order(KappaInd[, 1]), ]]
  colnames(Delta) <- paste0("Delta", 1:K)
  colnames(Tau) <- paste0("Tau", 1:K)
  UpsilonInd <- which(lower.tri(apply(matrix(1:K, ncol = 1), 1, function(x) paste0(paste0("Upsilon", 1:K, "_"), x)), diag = TRUE), arr.ind = TRUE)
  if (K == 1) colnames(Upsilon) <- apply(matrix(1:K, ncol = 1), 1, function(x) paste0(paste0("Upsilon", 1:K, "_"), x))[UpsilonInd[order(UpsilonInd[, 1]), ]][1]
  if (K > 1) colnames(Upsilon) <- apply(matrix(1:K, ncol = 1), 1, function(x) paste0(paste0("Upsilon", 1:K, "_"), x))[UpsilonInd[order(UpsilonInd[, 1]), ]]
  colnames(Psi) <- "Psi"
  colnames(Xi) <- paste0("Xi_", LambdaInd[, 3], "_", LambdaInd[, 2], "_", LambdaInd[, 1])
  colnames(Rho) <- "Rho"
  if (P == 0) Beta <- NULL
  if (P > 0) colnames(Beta) <- paste0("Beta", 1:P)
  if (CL == 0) Tau <- Delta <- Xi <- NULL
  Out <- list(Lambda = Lambda, Eta = Eta, Sigma2 = Sigma2, Kappa = Kappa, Delta = Delta, Tau = Tau, Upsilon = Upsilon, Psi = Psi, Xi = Xi, Rho = Rho, Beta = Beta)
  return(Out)
}



###Function for creating a data object that contains objects needed for ModelFit-----------------------------------
OutputDatObj <- function(DatObj) {

  ###Collect needed objects
  # DatObjOut <- list(M = DatObj$M,
  #                   Nu = DatObj$Nu,
  #                   AdjacentEdgesBoolean = DatObj$AdjacentEdgesBoolean,
  #                   W = DatObj$W,
  #                   EyeM = DatObj$EyeM,
  #                   EyeNu = DatObj$EyeNu,
  #                   OneM = DatObj$OneM,
  #                   OneN = DatObj$OneN,
  #                   OneNu = DatObj$OneNu,
  #                   YStarWide = DatObj$YStarWide,
  #                   Rho = DatObj$Rho,
  #                   FamilyInd = DatObj$FamilyInd,
  #                   ScaleY = DatObj$ScaleY,
  #                   YObserved = DatObj$YObserved,
  #                   ScaleDM = DatObj$ScaleDM,
  #                   Time = DatObj$Time,
  #                   TimeVec = DatObj$TimeVec,
  #                   YObserved = DatObj$YObserved,
  #                   tNu = DatObj$tNu,
  #                   t1 = DatObj$t1,
  #                   XThetaInd = DatObj$XThetaInd,
  #                   N = DatObj$N,
  #                   EyeN = DatObj$EyeN)
  DatObjOut <- DatObj
  return(DatObjOut)

}



###Function for creating a data augmentation object that contains objects needed for ModelFit----------------------
OutputDatAug <- function(DatAug) {

  ###Collect needed objects
  # DatAugOut <- list(NBelow = DatAug$NBelow,
  #                   NBelowList = DatAug$NBelowList,
  #                   TobitBooleanMat = DatAug$TobitBooleanMat,
  #                   YStarNonZeroList = DatAug$YStarNonZero)
  DatAugOut <- DatAug
  return(DatAugOut)

}



###Function for summarizing Metropolis objects post sampler--------------------------------------------------------
SummarizeMetropolis <- function(DatObj, MetrObj, MetropRcpp, McmcObj) {

  ###Set data object
  SpCorInd <- DatObj$SpCorInd
  IS <- DatObj$IS

  ###Set MCMC object
  NSims <- McmcObj$NSims

  ###Set Metropolis objects
  MetropPsi <- MetropRcpp$MetropPsi
  AcceptancePsi <- MetropRcpp$AcceptancePsi
  OriginalTuners <- MetrObj$OriginalTuners[1]
  AcceptancePct <- AcceptancePsi / NSims
  MetrSummary <- cbind(AcceptancePct, MetropPsi, OriginalTuners)
  rownames(MetrSummary) <- "Psi"
  colnames(MetrSummary) <- c("Acceptance", "PilotAdaptedTuners", "OriginalTuners")
  
  ###Add Rho
  if (SpCorInd == 0 & IS == 1) {
    MetropRho <- MetropRcpp$MetropRho
    AcceptanceRho <- MetropRcpp$AcceptanceRho
    OriginalTuners <- MetrObj$OriginalTuners
    AcceptancePct <- c(AcceptancePsi, AcceptanceRho) / NSims
    MetrSummary <- cbind(AcceptancePct, c(MetropPsi, MetropRho), OriginalTuners)
    rownames(MetrSummary) <- c("Psi", "Rho")
    colnames(MetrSummary) <- c("Acceptance", "PilotAdaptedTuners", "OriginalTuners")
  }

  ###Summarize and output
  return(MetrSummary)

}



###Verify the class of our regression object------------------------------------------------------------------------
#' is.spBFA
#'
#' \code{is.spBFA} is a general test of an object being interpretable as a
#' \code{\link{spBFA}} object.
#'
#' @param x object to be tested.
#'
#' @details The \code{\link{spBFA}} class is defined as the regression object that
#'  results from the \code{\link{spBFA}} regression function.
#'
#' @examples 
#' ###Load pre-computed results
#' data(reg.bfa_sp)
#' 
#' ###Test function
#' is.spBFA(reg.bfa_sp)
#'
#' @export
is.spBFA <- function(x) {
  identical(attributes(x)$class, "spBFA")
}