CheckInputs <- function(Y, DM, W, Time, Starting, Hypers, Tuning, MCMC, Family, Distance, Weights, Rho, ScaleY, ScaleDM) {

  ###Data dimensions
  N <- length(Y)
  M <- length(DM)
  Nu <- length(Time)

  ###Family
  if (!Family %in% c("normal", "probit", "tobit")) stop('Family: must be one of "normal", "probit" or "tobit"')

  ###Distance
  if (!Distance %in% c("euclidean", "circumference")) stop('Distance: must be one of "euclidean" or "circumference"')

  ###Weights
  if (!Weights %in% c("continuous", "binary")) stop('Weights: must be one of "continuous" or "binary"')

  ###Rho
  if (missing(Rho)) stop("Rho: missing")
  if (!is.scalar(Rho)) stop('Rho must be a scalar')
  if (is.na(Rho)) stop('Rho cannot be NA')
  if (!is.finite(Rho)) stop('Rho cannot be infinite')
  if (!((Rho < 1) & (Rho > 0))) stop('Rho must be in (0, 1)')

  ###ScaleY
  if (missing(ScaleY)) stop("ScaleY: missing")
  if (!is.scalar(ScaleY)) stop('ScaleY must be a scalar')
  if (is.na(ScaleY)) stop('ScaleY cannot be NA')
  if (!is.finite(ScaleY)) stop('ScaleY cannot be infinite')
  if (!(ScaleY > 0)) stop('ScaleY must be positive')

  ###ScaleDM
  if (missing(ScaleDM)) stop("ScaleDM: missing")
  if (!is.scalar(ScaleDM)) stop('ScaleDM must be a scalar')
  if (is.na(ScaleDM)) stop('ScaleDM cannot be NA')
  if (!is.finite(ScaleDM)) stop('ScaleDM cannot be infinite')
  if (!(ScaleDM > 0)) stop('ScaleDM must be positive')

  ###Data checks for Y
  if (!is.numeric(Y)) stop('Y must be a vector')
  if (length(Y) != N) stop(paste0('Y must have length ', N))
  if (any(is.na(Y))) stop("Y may have no missing values")
  if (any(!is.finite(Y))) stop("Y must have strictly finite entries")
  if ((Family == "probit") & ((sum(Y == 1) + sum(Y == 0)) != N)) stop('Y: for "probit" observed data must be binary')
  if ((Family == "tobit") & (any(Y < 0))) stop('Y: for "tobit" observed data must be non-negative')

  ###Data checks for DM
  if (!is.numeric(DM)) stop('DM must be a vector')
  if (length(DM) != M) stop(paste0('DM must have length ', M))
  if (any(is.na(DM))) stop("DM may have no missing values")
  if (any(!is.finite(DM))) stop("DM must have strictly finite entries")

  ###Data checks for W
  if (!is.matrix(W)) stop('W must be a matrix')
  if (!dim(W)[1] == M) stop(paste0('W must be a ',M ,' x ', M, ' dimensional matrix'))
  if (!dim(W)[2] == M) stop('W must be square')
  if (sum(!((W) == t(W))) > 0) stop('W must be symmetric')
  if (length(table(W)) > 2) stop('W must only contain binaries (i.e. 0\'s or 1\'s)')
  if (any(diag(W) != 0)) stop('W must have only zeros on the diagonal')

  ###Data checks for Time
  if (!is.numeric(Time)) stop('Time must be a vector')
  if (length(Time) != Nu) stop(paste0('Time must have length ', Nu))
  if (any(is.na(Time))) stop("Time may have no missing values")
  if (any(!is.finite(Time))) stop("Time must have strictly finite entries")
  if (is.unsorted(Time)) stop('Time vector is not in increasing order')
  if (!all(Time >= 0)) stop('Time vector has at least one negative point')

  ###Verify dimensions
  M_W <- dim(W)[1]
  if (M != M_W) stop('DM and W have contradictory dimensions')
  if ((N / M) != Nu) stop('Time, DM and Y have contradictory dimensions')

  ###Hypers
  if (!is.null(Hypers)) {
    if (!is.list(Hypers)) stop('Hypers must be a list')
    if (!all(names(Hypers) %in% c("Delta", "Sigma", "Alpha"))) stop('Hypers: Can only contain lists with names "Delta", "Sigma" and "Alpha"')

    ###If delta hyperparameters are provided
    if ("Delta" %in% names(Hypers)) {
      if (!is.list(Hypers$Delta)) stop('Hypers: "Delta" must be a list')
      if (!"Kappa2" %in% names(Hypers$Delta)) stop('Hypers: "Kappa2" value missing')
      if (!is.scalar(Hypers$Delta$Kappa2)) stop('Hypers: "Kappa2" must be a scalar')
      if (is.na(Hypers$Delta$Kappa2)) stop('Hypers: "Kappa2" cannot be NA')
      if (!is.finite(Hypers$Delta$Kappa2)) stop('Hypers: "Kappa2" cannot be infinite')
      if (Hypers$Delta$Kappa2 <= 0) stop('Hypers: "Kappa2" must be strictly positive')
    }

    ###If Sigma hyperparameters are provided
    if ("Sigma" %in% names(Hypers)) {
      if (!is.list(Hypers$Sigma)) stop('Hypers: "Sigma" must be a list')
      if (!"Xi" %in% names(Hypers$Sigma)) stop('Hypers: "Xi" value missing')
      if (!is.scalar(Hypers$Sigma$Xi)) stop('Hypers: "Xi" must be a scalar')
      if (is.na(Hypers$Sigma$Xi)) stop('Hypers: "Xi" cannot be NA')
      if (!is.finite(Hypers$Sigma$Xi)) stop('Hypers: "Xi" cannot be infinite')
      if (Hypers$Sigma$Xi < 3) stop('Hypers: "Xi" must be greater than or equal to 5')
      if (!"Psi" %in% names(Hypers$Sigma)) stop('Hypers: "Psi" value missing')
      if (!is.matrix(Hypers$Sigma$Psi)) stop('Hypers: "Psi" must be a matrix')
      if (!dim(Hypers$Sigma$Psi)[1] == 5) stop('Hypers: "Psi" must be 5 dimensional')
      if (!all(!is.na(Hypers$Sigma$Psi))) stop('Hypers: "Psi" cannot have missing values')
      if (!all(is.finite(Hypers$Sigma$Psi))) stop('Hypers: "Psi" cannot have infinite values')
      if (!dim(Hypers$Sigma$Psi)[2] == 5) stop('Hypers: "Psi" must be square')
      if (sum( !( (Hypers$Sigma$Psi) == t(Hypers$Sigma$Psi) ) ) > 0) stop('Hypers: "Psi" must be symmetric')
      if ((det(Hypers$Sigma$Psi) - 0) < 0.00001) stop('Hypers: "Psi" is close to singular')
    }

    ###If Alpha hyperparameters are provided
    if ("Alpha" %in% names(Hypers)) {
      if (!is.list(Hypers$Alpha)) stop('Hypers: "Alpha" must be a list')
      if (!"AAlpha" %in% names(Hypers$Alpha)) stop('Hypers: "AAlpha" value missing')
      if (!is.scalar(Hypers$Alpha$AAlpha)) stop('Hypers: "AAlpha" must be a scalar')
      if (is.na(Hypers$Alpha$AAlpha)) stop('Hypers: "AAlpha" cannot be NA')
      if (!is.finite(Hypers$Alpha$AAlpha)) stop('Hypers: "AAlpha" cannot be infinite')
      if (!"BAlpha" %in% names(Hypers$Alpha)) stop('Hypers: "BAlpha" value missing')
      if (!is.scalar(Hypers$Alpha$BAlpha)) stop('Hypers: "BAlpha" must be a scalar')
      if (is.na(Hypers$Alpha$BAlpha)) stop('Hypers: "BAlpha" cannot be NA')
      if (!is.finite(Hypers$Alpha$BAlpha)) stop('Hypers: "BAlpha" cannot be infinite')
      if (Hypers$Alpha$AAlpha < 0) stop('Hypers: "AAlpha" must be non-negative')
      if (Hypers$Alpha$BAlpha <= 0) stop('Hypers: "BAlpha" must be strictly positive')
      if (Hypers$Alpha$BAlpha < Hypers$Alpha$AAlpha) stop('Hypers: "BAlpha" must be greater than "AAlpha"')
    }

  ###End Hyperparameters
  }

  ###Starting Values
  if (!is.null(Starting)) {
    if (!is.list(Starting)) stop('Starting must be a list')
    if (!all(names(Starting) %in% c("Delta", "Sigma", "Alpha"))) stop('Starting: Can only contain objects with names "Delta", "Sigma" and "Alpha"')

    ###If delta starting values is provided
    if ("Delta" %in% names(Starting)) {
      if (!is.numeric(Starting$Delta)) stop('Starting: "Delta" must be a vector')
      if (length(Starting$Delta) != 5) stop('Starting: "Delta" must be length 5')
      if (!all(!is.na(Starting$Delta))) stop('Starting: "Delta" cannot have missing values')
      if (!all(is.finite(Starting$Delta))) stop('Starting: "Delta" cannot have infinite values')
    }

    ###If Sigma starting values is provided
    if ("Sigma" %in% names(Starting)) {
      if (!is.matrix(Starting$Sigma)) stop('Starting: "Sigma" must be a matrix')
      if (!dim(Starting$Sigma)[1] == 5) stop('Starting: "Sigma" must be 5 dimensional')
      if (!dim(Starting$Sigma)[2] == 5) stop('Starting: "Sigma" must be square')
      if (!all(!is.na(Starting$Sigma))) stop('Starting: "Sigma" cannot have missing values')
      if (!all(is.finite(Starting$Sigma))) stop('Starting: "Sigma" cannot have infinite values')
      if (sum( !( (Starting$Sigma) == t(Starting$Sigma) ) ) > 0) stop('Starting: "Sigma" must be symmetric')
      if ((det(Starting$Sigma) - 0) < 0.0000000001) stop('Starting: "Sigma" is close to singular')
    }

    ###If Alpha starting values is provided
    if ("Alpha" %in% names(Starting)) {
      if (!is.scalar(Starting$Alpha)) stop('Starting: "Alpha" must be a scalar')
      if (is.na(Starting$Alpha)) stop('Starting: "Alpha" cannot be NA')
      if (!is.finite(Starting$Alpha)) stop('Starting: "Alpha" cannot be infinite')
      # I make sure that Alpha is in (AAlpha, BAlpha) in CreatePara();
    }

  ###End Starting Values
  }

  ###Tuning Values
  if (!is.null(Tuning)) {
    if (!is.list(Tuning)) stop('Tuning must be a list')
    if (!all(names(Tuning) %in% c("Lambda0Vec", "Lambda1Vec", "EtaVec", "Alpha"))) stop('Tuning: Can only contain objects with names "Lambda0Vec", "Lambda1Vec", "EtaVec" or "Alpha"')

    ###If Lambda0Vec tuning values are provided
    if ("Lambda0Vec" %in% names(Tuning)) {
      if (!is.numeric(Tuning$Lambda0Vec)) stop('Tuning: "Lambda0Vec" must be a vector')
      if (length(Tuning$Lambda0Vec) != M) stop(paste0('Tuning: "Lambda0Vec" must have length ', M))
      if (!all(!is.na(Tuning$Lambda0Vec))) stop('Tuning: "Lambda0Vec" cannot have missing values')
      if (!all(is.finite(Tuning$Lambda0Vec))) stop('Tuning: "Lambda0Vec" cannot have infinite values')
      if (any(Tuning$Lambda0Vec < 0)) stop('Tuning: "Lambda0Vec" must have non-negative components')
    }

    ###If Lambda1Vec tuning values are provided
    if ("Lambda1Vec" %in% names(Tuning)) {
      if (!is.numeric(Tuning$Lambda1Vec)) stop('Tuning: "Lambda1Vec" must be a vector')
      if (length(Tuning$Lambda1Vec) != M) stop(paste0('Tuning: "Lambda1Vec" must have length ', M))
      if (!all(!is.na(Tuning$Lambda1Vec))) stop('Tuning: "Lambda1Vec" cannot have missing values')
      if (!all(is.finite(Tuning$Lambda1Vec))) stop('Tuning: "Lambda1Vec" cannot have infinite values')
      if (any(Tuning$Lambda1Vec < 0)) stop('Tuning: "Lambda1Vec" must have non-negative components')
    }

    ###If EtaVec tuning values are provided
    if ("EtaVec" %in% names(Tuning)) {
      if (!is.numeric(Tuning$EtaVec)) stop('Tuning: "EtaVec" must be a vector')
      if (length(Tuning$EtaVec) != M) stop(paste0('Tuning: "EtaVec" must have length ', M))
      if (!all(!is.na(Tuning$EtaVec))) stop('Tuning: "EtaVec" cannot have missing values')
      if (!all(is.finite(Tuning$EtaVec))) stop('Tuning: "EtaVec" cannot have infinite values')
      if (any(Tuning$EtaVec < 0)) stop('Tuning: "EtaVec" must have non-negative components')
    }

    ###If Alpha tuning value is provided
    if ("Alpha" %in% names(Tuning)) {
      if (!is.scalar(Tuning$Alpha)) stop('Tuning: "Alpha" must be a scalar')
      if (is.na(Tuning$Alpha)) stop('Tuning: "Alpha" cannot be NA')
      if (!is.finite(Tuning$Alpha)) stop('Tuning: "Alpha" cannot be infinite')
      if (Tuning$Alpha < 0) stop('Tuning: "Alpha" must be non-negative')
    }

  ###End Tuning Values
  }

  ###MCMC Values
  if (!is.null(MCMC)) {
    if (!is.list(MCMC)) stop('MCMC must be a list')
    if (!all(names(MCMC) %in% c("NBurn", "NSims", "NThin", "NPilot"))) stop('MCMC: Can only contain objects with names "NBurn", "NSims", "NThin" and "NPilot"')

    ###If NBurn is provided
    if ("NBurn" %in% names(MCMC)) {
      if (!is.scalar(MCMC$NBurn)) stop('MCMC: "NBurn" must be a scalar')
      if (is.na(MCMC$NBurn)) stop('MCMC: "NBurn" cannot be NA')
      if (!is.finite(MCMC$NBurn)) stop('MCMC: "NBurn" cannot be infinite')
      if (!is.wholenumber(MCMC$NBurn) | MCMC$NBurn < 0) stop('MCMC: "NBurn" must be a non-negative integer')
      if (MCMC$NBurn < 100) stop('MCMC: "NBurn" must be at least 100')
    }

    ###If NSims is provided
    if ("NSims" %in% names(MCMC)) {
      if (!is.scalar(MCMC$NSims)) stop('MCMC: "NSims" must be a scalar')
      if (is.na(MCMC$NSims)) stop('MCMC: "NSims" cannot be NA')
      if (!is.finite(MCMC$NSims)) stop('MCMC: "NSims" cannot be infinite')
      if (!is.wholenumber(MCMC$NSims) | MCMC$NSims <= 0) stop('MCMC: "NSims" must be a positive integer')
      if (MCMC$NSims < 100) stop('MCMC: "NSims" must be at least 100')
    }

    ###If NThin is provided
    if ("NThin" %in% names(MCMC)) {
      if (!is.scalar(MCMC$NThin)) stop('MCMC: "NThin" must be a scalar')
      if (is.na(MCMC$NThin)) stop('MCMC: "NThin" cannot be NA')
      if (!is.finite(MCMC$NThin)) stop('MCMC: "NThin" cannot be infinite')
      if (!is.wholenumber(MCMC$NThin) | MCMC$NThin <= 0) stop('MCMC: "NThin" must be a positive integer')
      # if (!is.wholenumber(MCMC$NSims / MCMC$NThin)) stop('MCMC: "NThin" must be a factor of "NSims"') enforced in CreateMCMC();
    }

    ###If NPilot is provided
    if ("NPilot" %in% names(MCMC)) {
      if (!is.scalar(MCMC$NPilot)) stop('MCMC: "NPilot" must be a scalar')
      if (is.na(MCMC$NPilot)) stop('MCMC: "NPilot" cannot be NA')
      if (!is.finite(MCMC$NPilot)) stop('MCMC: "NPilot" cannot be infinite')
      if (!is.wholenumber(MCMC$NPilot) | MCMC$NPilot < 0) stop('MCMC: "NPilot" must be a positive integer')
      # if (!is.wholenumber(MCMC$NBurn / MCMC$NPilot)) stop('MCMC: "NPilot" must be a factor of "NBurn"') enforced in CreateMCMC();
    }

  ###End MCMC Values
  }

}

###Helper Functions
is.scalar <- function(x) ((is.numeric(x)) & (length(x) == 1))
is.wholenumber <- function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol
