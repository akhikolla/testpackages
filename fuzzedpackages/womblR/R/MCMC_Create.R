###Function for reading in sampler inputs and creating a list object that contains all relavent data objects--------------------
CreateDatObj <- function(Y, DM, W, Time, Rho, ScaleY, ScaleDM, TemporalStructure, Family, Distance, Weights) {

  ###Data objects
  YObserved <- Y / ScaleY #scale observed data
  N <- length(YObserved)  #total observations
  M <- dim(W)[1] #number of spatial locations
  Nu <- N / M #number of visits
  NTheta <- 3 * Nu #number of theta parameters
  YObserved <- matrix(YObserved, ncol = 1)
  DM <- matrix(DM, ncol = 1)

  ###Temporal distance matrix
  TimeDist <- abs( outer(Time, Time, "-"))

  ###Dynamic data object (updated with data augmentation)
  YStarWide <- matrix(YObserved, nrow = M, ncol = Nu)

  ###Matrix Objects
  OneM <- matrix(rep(1, M), nrow = M , ncol = 1)
  OneNu <- matrix(rep(1, Nu), nrow = Nu, ncol = 1)
  EyeM <- diag(M)
  EyeNu <- diag(Nu)
  EyeNTheta <- diag(3 * Nu)
  Eye3 <- diag(3)
  ZDelta <- kronecker(OneNu, Eye3)

  ###Create Z vector
  if (Distance == "euclidean") Dist <- function(x, y) abs(x - y)
  if (Distance == "circumference") Dist <- function(x, y) pmin(abs(x - y), (360 - pmax(x, y) + pmin(x, y))) #circumference of optic nerve
  DM_Grid <- expand.grid(DM, DM)
  Z_Vector <- Dist(DM_Grid[ , 1], DM_Grid[ , 2])
  Z_Matrix <- matrix(Z_Vector, nrow = dim(W)[1], ncol = dim(W)[1], byrow = TRUE)
  AdjacentEdgesBoolean <- (W == 1) & (!lower.tri(W))
  Z <- matrix(Z_Matrix[AdjacentEdgesBoolean] / ScaleDM, ncol = 1)
  AdjacentEdgesBoolean <- matrix(which(AdjacentEdgesBoolean) - 1, ncol = 1)

  ###Assign temporal correlation structure
  if (TemporalStructure == "exponential") TempCorInd <- 0
  if (TemporalStructure == "ar1") TempCorInd <- 1

  ###Family indicator
  if (Family == "normal") FamilyInd <- 0
  if (Family == "probit") FamilyInd <- 1
  if (Family == "tobit") FamilyInd <- 2

  ###Weights indicator
  if (Weights == "continuous") WeightsInd <- 0
  if (Weights == "binary") WeightsInd <- 1

  ###Make parameters global
  DatObj <- list()
  DatObj$YObserved <- YObserved
  DatObj$ScaleY <- ScaleY
  DatObj$ScaleDM <- ScaleDM
  DatObj$YStarWide <- YStarWide
  DatObj$DM <- DM
  DatObj$W <- W
  DatObj$TimeDist <- TimeDist
  DatObj$Rho <- Rho
  DatObj$N <- N
  DatObj$M <- M
  DatObj$Nu <- Nu
  DatObj$NTheta <- NTheta
  DatObj$OneM <- OneM
  DatObj$OneNu <- OneNu
  DatObj$EyeM <- EyeM
  DatObj$EyeNu <- EyeNu
  DatObj$EyeNTheta <- EyeNTheta
  DatObj$Eye3 <- Eye3
  DatObj$ZDelta <- ZDelta
  DatObj$Z <- Z
  DatObj$AdjacentEdgesBoolean <- AdjacentEdgesBoolean
  DatObj$TempCorInd <- TempCorInd
  DatObj$FamilyInd <- FamilyInd
  DatObj$WeightsInd <- WeightsInd
  DatObj$Time <- Time
  return(DatObj)

}



###Function to create Hyperparameter Object------------------------------------------------------------------------------------
CreateHyPara <- function(Hypers, DatObj) {

  ###Set data objects
  TimeDist <- DatObj$TimeDist
  TempCorInd <- DatObj$TempCorInd

  ###Which parameters are user defined?
  UserHypers <- names(Hypers)

  ###Set hyperparameters for Delta
  if ("Delta" %in% UserHypers) {
    MuDelta <- matrix(Hypers$Delta$MuDelta, nrow = 3)
    OmegaDeltaInv <- solve(Hypers$Delta$OmegaDelta)
    OmegaDeltaInvMuDelta <- OmegaDeltaInv %*% MuDelta
  }
  if (!"Delta" %in% UserHypers) {
    MuDelta <- matrix(c(3, 0, 0), nrow = 3)
    OmegaDeltaInv <- diag(c(0.001, 0.001, 1))
    OmegaDeltaInvMuDelta <- OmegaDeltaInv %*% MuDelta
  }

  ###Set hyperparameters for T
  if ("T" %in% UserHypers) {
    Xi <- Hypers$T$Xi
    Psi <- Hypers$T$Psi
  }
  if (!"T" %in% UserHypers) {
    Xi <- 3 + 1
    Psi <- diag(3)
  }

  ###Set hyperparameters for Phi
  if ("Phi" %in% UserHypers) {
    APhi <- Hypers$Phi$APhi
    BPhi <- Hypers$Phi$BPhi
  }
  if (!"Phi" %in% UserHypers) {
    minDiff <- min( TimeDist[ TimeDist > 0 ] )
    maxDiff <- max( TimeDist[ TimeDist > 0 ] )
    if (TempCorInd == 0) { # exponential
      BPhi <- -log(0.01) / minDiff #shortest diff goes down to 1%
      APhi <- -log(0.95) / maxDiff #longest diff goes up to 95%
    }
    if (TempCorInd == 1) { # ar1
      APhi <- 0.01 ^ (1 / minDiff) #shortest diff goes down to 1%
      BPhi <- 0.95 ^ (1 / maxDiff) #longest diff goes up to 95%
    }
  }

  ###Create object for hyperparameters
  HyPara <- list()
  HyPara$OmegaDeltaInvMuDelta <- OmegaDeltaInvMuDelta
  HyPara$OmegaDeltaInv <- OmegaDeltaInv
  HyPara$Xi <- Xi
  HyPara$Psi <- Psi
  HyPara$APhi <- APhi
  HyPara$BPhi <- BPhi
  return(HyPara)

}



###Function for creating an object containing relevant Metropolis information---------------------------------------------------
CreateMetrObj <- function(Tuning, DatObj) {

  ###Set Data Objects
  Nu <- DatObj$Nu

  ###Which parameters are user defined?
  UserTuners <- names(Tuning)

  ###Set tuning parameters for Theta2
  if ("Theta2" %in% UserTuners) MetropTheta2 <- Tuning$Theta2
  if (!"Theta2" %in% UserTuners) MetropTheta2 <- rep(1, Nu)

  ###Set tuning parameters for Theta3
  if ("Theta3" %in% UserTuners) MetropTheta3 <- Tuning$Theta3
  if (!"Theta3" %in% UserTuners) MetropTheta3 <- rep(1, Nu)

  ###Set tuning parameter for Phi
  if ("Phi" %in% UserTuners) MetropPhi <- Tuning$Phi
  if (!"Phi" %in% UserTuners) MetropPhi <- 1

  ###Set acceptance rate counters
  AcceptanceTheta2 <- AcceptanceTheta3 <- rep(0, Nu)
  AcceptancePhi <- 0

  ###Return metropolis object
  MetrObj <- list()
  MetrObj$MetropTheta2 <- MetropTheta2
  MetrObj$AcceptanceTheta2 <- AcceptanceTheta2
  MetrObj$MetropTheta3 <- MetropTheta3
  MetrObj$AcceptanceTheta3 <- AcceptanceTheta3
  MetrObj$MetropPhi <- MetropPhi
  MetrObj$AcceptancePhi <- AcceptancePhi
  return(MetrObj)

}



###Function for creating inital parameter object-------------------------------------------------------------------------------
CreatePara <- function(Starting, DatObj, HyPara) {

  ###Set data objects
  W <- DatObj$W
  Rho <- DatObj$Rho
  TimeDist <- DatObj$TimeDist
  Nu <- DatObj$Nu
  M <- DatObj$M
  ZDelta <- DatObj$ZDelta
  Z <- DatObj$Z
  AdjacentEdgesBoolean <- DatObj$AdjacentEdgesBoolean
  EyeM <- DatObj$EyeM
  TempCorInd <- DatObj$TempCorInd
  EyeNTheta <- DatObj$EyeNTheta
  OneNu <- DatObj$OneNu
  WeightsInd <- DatObj$WeightsInd

  ###Set hyperparameter objects
  APhi <- HyPara$APhi
  BPhi <- HyPara$BPhi

  ###Which parameters are user defined?
  UserStarters <- names(Starting)

  ###Set initial value of Delta
  if ("Delta" %in% UserStarters) Delta <- matrix(Starting$Delta, nrow = 3)
  if (!"Delta" %in% UserStarters) Delta <- matrix(c(3, 0, 0), nrow = 3, ncol = 1)

  ###Set intial value of T
  if ("T" %in% UserStarters) T <- Starting$T
  if (!"T" %in% UserStarters) T <- diag(3)

  ###Set initial values of Phi
  if ("Phi" %in% UserStarters) {
    Phi <- Starting$Phi
    if ((Phi <= APhi) | (Phi >= BPhi)) stop('Starting: "Phi" must be in the interval (APhi, BPhi)')
  }
  if (!"Phi" %in% UserStarters) Phi <- mean(c(APhi, BPhi))

  ###Set inital value of Theta (both matrix and vector form)
  VecTheta <- ZDelta %*% Delta
  Theta <- matrix(VecTheta, nrow = 3, ncol = Nu)

  ###Transform to level 1 parameters
  Mu <- Theta[1 , ]
  Tau2 <- exp( Theta[2 , ] ) ^ 2
  Alpha <- exp( Theta[ 3 , ] )

  ###Create covariance arrays that can be converted to arma::cubes
  WAlphas <- WAlphaCube(Alpha, Z, W, M, Nu, WeightsInd)
  JointCovariances <- JointCovarianceCube(WAlphas, Tau2, EyeM, Rho, M, Nu)
  RootiLikelihoods <- RootiLikelihoodCube(JointCovariances, EyeM, M, Nu)

  ###Prior covariance objects
  SIGMAPhi <- SIGMA(Phi, TempCorInd, TimeDist, Nu)
  SIGMAPhiInv <- CholInv(SIGMAPhi)
  TInv <- Inv3(T)
  CovTheta <- kronecker(SIGMAPhi, T)
  CovThetaInv <- kronecker(SIGMAPhiInv, TInv)
  RootiTheta <- GetRooti(CovTheta, EyeNTheta)
  MMat <- Delta %*% t(OneNu)
  MeanTheta <- ZDelta %*% Delta

  ###Save parameter objects
  Para <- list()
  Para$Mu <- Mu
  Para$Tau2 <- Tau2
  Para$Alpha <- Alpha
  Para$WAlphas <- WAlphas
  Para$JointCovariances <- JointCovariances
  Para$RootiLikelihoods <- RootiLikelihoods
  Para$VecTheta <- VecTheta
  Para$Theta <- Theta
  Para$Delta <- Delta
  Para$MeanTheta <- MeanTheta
  Para$T <- T
  Para$TInv <- TInv
  Para$Phi <- Phi
  Para$SIGMAPhi <- SIGMAPhi
  Para$SIGMAPhiInv <- SIGMAPhiInv
  Para$CovThetaInv <- CovThetaInv
  Para$RootiTheta <- RootiTheta
  Para$MMat <- MMat
  return(Para)

}



###Function that creates the data augmentation (i.e. Tobit) booleans------------------------------------------------------------
CreateDatAug <- function(DatObj) {

  ###Set data object
  YObserved <- DatObj$YObserved
  FamilyInd <- DatObj$FamilyInd
  YStarWide <- DatObj$YStarWide
  M <- DatObj$M
  Nu <- DatObj$Nu

  ###Initialize Data Augmentation Object
  DatAug <- NULL

  ###Normal objects
  if (FamilyInd == 0) {
    DatAug$NBelow <- 0
    DatAug$NAbove <- 0
    DatAug$TobitIndeces <- matrix(c(0,0,0,0),2,2)
    DatAug$ProbitIndeces <- matrix(c(0,0,0,0),2,2)
  }

  ###Probit objects
  if (FamilyInd == 1) {
    TobitBoolean <- YObserved <= 0
    WhichBelow <- which(TobitBoolean)
    NBelow <- length(WhichBelow)
    TobitBooleanMat <- matrix(TobitBoolean, nrow = M, ncol = Nu)
    YStarBelow <- list()
    for (i in 1 : Nu) YStarBelow[[i]] <- YStarWide[!TobitBooleanMat[,i], i]
    NBelowList <- unlist(lapply(YStarBelow, f<-function(x) M - length(x)))
    TobitIndeces <- which(TobitBooleanMat, arr.ind = TRUE)
    TobitIndeces <- TobitIndeces - 1
    ProbitBoolean <- YObserved > 0
    WhichAbove <- which(ProbitBoolean)
    NAbove <- length(WhichAbove)
    ProbitBooleanMat <- matrix(ProbitBoolean, nrow = M, ncol = Nu)
    YStarAbove <- list()
    for (i in 1 : Nu) YStarAbove[[i]] <- YStarWide[!ProbitBooleanMat[,i], i]
    NAboveList <- unlist(lapply(YStarAbove, f<-function(x) M - length(x)))
    ProbitIndeces <- which(ProbitBooleanMat, arr.ind = TRUE)
    ProbitIndeces <- ProbitIndeces - 1

    ###Save objects
    DatAug <- list()
    DatAug$NBelow <- NBelow
    DatAug$NAbove <- NAbove
    DatAug$TobitIndeces <- TobitIndeces
    DatAug$ProbitIndeces <- ProbitIndeces
  }

  ###Tobit objects
  if (FamilyInd == 2) {
    TobitBoolean <- YObserved <= 0
    WhichBelow <- which(TobitBoolean)
    NBelow <- length(WhichBelow)
    TobitBooleanMat <- matrix(TobitBoolean, nrow = M, ncol = Nu)
    YStarNonZero <- list()
    for (i in 1 : Nu) YStarNonZero[[i]] <- YStarWide[!TobitBooleanMat[,i], i]
    NBelowCount <- unlist(lapply(YStarNonZero, f<-function(x) M - length(x)))
    TobitIndeces <- which(TobitBooleanMat, arr.ind = TRUE)
    TobitIndeces <- TobitIndeces - 1
    # ZDatAug <- model.matrix(~-1 + as.factor(TobitIndeces[,2]))
    # attributes(ZDatAug) <- NULL
    # ZDatAug <- structure(ZDatAug, class = "matrix", dim = c(NBelow, Nu))
    # WDatAug <- array(FALSE, dim = c(M, M, Nu))
    # for (i in 1:NBelow) {
    #   Visit <- TobitIndeces[i, 2] + 1
    #   Location <- TobitIndeces[i, 1] + 1
    #   WDatAug[ , Location, Visit] <- rep(TRUE, M)
    # }
    # WDatAug <- matrix(which(WDatAug) - 1, ncol = 1)

    ###Save objects
    DatAug <- list()
    DatAug$WhichBelow <- WhichBelow
    DatAug$NBelow <- NBelow
    DatAug$TobitBooleanMat <- TobitBooleanMat
    DatAug$YStarNonZero <- YStarNonZero
    DatAug$NBelowCount <- NBelowCount
    DatAug$TobitIndeces <- TobitIndeces
    DatAug$ProbitIndeces <- matrix(c(0,0,0,0),2,2)
    DatAug$NAbove <- 2
    # DatAug$ZDatAug <- ZDatAug
    # DatAug$WDatAug <- WDatAug
  }

  return(DatAug)

}



###Function that creates inputs for MCMC sampler--------------------------------------------------------------------------------
CreateMcmc <- function(MCMC, DatObj) {

  ###Set data objects
  Nu <- DatObj$Nu

  ###Which parameters are user defined?
  UserMCMC <- names(MCMC)

  ###Set MCMC objects
  if ("NBurn" %in% UserMCMC) NBurn <- MCMC$NBurn
  if (!"NBurn" %in% UserMCMC) NBurn <- 10000
  if ("NSims" %in% UserMCMC) NSims <- MCMC$NSims
  if (!"NSims" %in% UserMCMC) NSims <- 100000
  if ("NThin" %in% UserMCMC) NThin <- MCMC$NThin
  if (!"NThin" %in% UserMCMC) NThin <- 10
  if ("NPilot" %in% UserMCMC) NPilot <- MCMC$NPilot
  if (!"NPilot" %in% UserMCMC) NPilot <- 20

  ###One last check of MCMC user inputs
  is.wholenumber <- function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol
  if (!is.wholenumber(NSims / NThin)) stop('MCMC: "NThin" must be a factor of "NSims"')
  if (!is.wholenumber(NBurn / NPilot)) stop('MCMC: "NPilot" must be a factor of "NBurn"')

  ###Create MCMC objects
  NTotal <- NBurn + NSims
  WhichKeep <- NBurn + ( 1 : (NSims / NThin) ) * NThin
  NKeep <- length(WhichKeep)

  ###Pilot adaptation objects
  WhichPilotAdapt <- ( 1 : NPilot ) * NBurn / NPilot
  PilotAdaptDenominator <- WhichPilotAdapt[1]

  ###Burn-in progres bar
  BarLength <- 50 #Burn-in bar length (arbitrary)
  BurnInProgress <- seq(1 / BarLength, 1, 1 / BarLength)
  WhichBurnInProgress <- sapply(BurnInProgress, function(x) tail(which(1 : NBurn <= x * NBurn), 1))

  ###Progress output objects
  SamplerProgress <- seq(0.1, 1.0, 0.1) #Intervals of progress update (arbitrary)
  WhichSamplerProgress <- sapply(SamplerProgress, function(x) tail(which(1 : NSims <= x * NSims), 1)) + NBurn
  WhichBurnInProgressInt <- sapply(SamplerProgress, function(x) tail(which(1 : NBurn <= x * NBurn), 1))


  ###Save objects
  MCMC <- list()
  MCMC$NBurn <- NBurn
  MCMC$NSims <- NSims
  MCMC$NThin <- NThin
  MCMC$NPilot <- NPilot
  MCMC$NTotal <- NTotal
  MCMC$WhichKeep <- WhichKeep
  MCMC$NKeep <- NKeep
  MCMC$WhichPilotAdapt <- WhichPilotAdapt
  MCMC$PilotAdaptDenominator <- PilotAdaptDenominator
  MCMC$BurnInProgress <- BurnInProgress
  MCMC$WhichBurnInProgress <- WhichBurnInProgress
  MCMC$WhichBurnInProgressInt <- WhichBurnInProgressInt
  MCMC$BarLength <- BarLength
  MCMC$WhichSamplerProgress <- WhichSamplerProgress
  return(MCMC)

}



###Function that creates a storage object for raw samples-----------------------------------------------------------------------
CreateStorage <- function(DatObj, McmcObj) {

  ###Set data objects
  Nu <- DatObj$Nu

  ###Set MCMC objects
  NKeep <- McmcObj$NKeep

  ###Create storage object
  Out <- matrix(nrow = (3 * Nu + 3 + 6 + 1), ncol = NKeep)
  return(Out)

}

