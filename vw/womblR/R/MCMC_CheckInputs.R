CheckInputs <- function(Y, DM, W, Time, Starting, Hypers, Tuning, MCMC, Family, TemporalStructure, Distance, Weights, Rho, ScaleY, ScaleDM) {

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

  ###Temporal structure
  if (!TemporalStructure %in% c("exponential", "ar1")) stop('TemporalStructure: must be one of "exponential" or "ar1"')

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
    if (!all(names(Hypers) %in% c("Delta", "T", "Phi"))) stop('Hypers: Can only contain lists with names "Delta", "T" and "Phi"')

    ###If delta hyperparameters are provided
    if ("Delta" %in% names(Hypers)) {
      if (!is.list(Hypers$Delta)) stop('Hypers: "Delta" must be a list')
      if (!"MuDelta" %in% names(Hypers$Delta)) stop('Hypers: "MuDelta" value missing')
      if (!is.numeric(Hypers$Delta$MuDelta)) stop('Hypers: "MuDelta" must be a vector')
      if (length(Hypers$Delta$MuDelta) != 3) stop('Hypers: "MuDelta" must be length 3')
      if (!all(!is.na(Hypers$Delta$MuDelta))) stop('Hypers: "MuDelta" cannot have missing values')
      if (!all(is.finite(Hypers$Delta$MuDelta))) stop('Hypers: "MuDelta" cannot have infinite values')
      if (!"OmegaDelta" %in% names(Hypers$Delta)) stop('Hypers: "OmegaDelta" value missing')
      if (!is.matrix(Hypers$Delta$OmegaDelta)) stop('Hypers: "OmegaDelta" must be a matrix')
      if (!dim(Hypers$Delta$OmegaDelta)[1] == 3) stop('Hypers: "OmegaDelta" must be 3 dimensional')
      if (!all(!is.na(Hypers$Delta$OmegaDelta))) stop('Hypers: "OmegaDelta" cannot have missing values')
      if (!all(is.finite(Hypers$Delta$OmegaDelta))) stop('Hypers: "OmegaDelta" cannot have infinite values')
      if (!dim(Hypers$Delta$OmegaDelta)[2] == 3) stop('Hypers: "OmegaDelta" must be square')
      if (sum( !( (Hypers$Delta$OmegaDelta) == t(Hypers$Delta$OmegaDelta) ) ) > 0) stop('Hypers: "OmegaDelta" must be symmetric')
      if ((det(Hypers$Delta$OmegaDelta) - 0) < 0.00001) stop('Hypers: "OmegaDelta" is close to singular')
      if (!all(!is.na(Hypers$Delta$OmegaDelta))) stop('Hypers: "OmegaDelta" cannot have missing values')
    }

    ###If T hyperparameters are provided
    if ("T" %in% names(Hypers)) {
      if (!is.list(Hypers$T)) stop('Hypers: "T" must be a list')
      if (!"Xi" %in% names(Hypers$T)) stop('Hypers: "Xi" value missing')
      if (!is.scalar(Hypers$T$Xi)) stop('Hypers: "Xi" must be a scalar')
      if (is.na(Hypers$T$Xi)) stop('Hypers: "Xi" cannot be NA')
      if (!is.finite(Hypers$T$Xi)) stop('Hypers: "Xi" cannot be infinite')
      if (Hypers$T$Xi < 3) stop('Hypers: "Xi" must be greater than or equal to 3')
      if (!"Psi" %in% names(Hypers$T)) stop('Hypers: "Psi" value missing')
      if (!is.matrix(Hypers$T$Psi)) stop('Hypers: "Psi" must be a matrix')
      if (!dim(Hypers$T$Psi)[1] == 3) stop('Hypers: "Psi" must be 3 dimensional')
      if (!all(!is.na(Hypers$T$Psi))) stop('Hypers: "Psi" cannot have missing values')
      if (!all(is.finite(Hypers$T$Psi))) stop('Hypers: "Psi" cannot have infinite values')
      if (!dim(Hypers$T$Psi)[2] == 3) stop('Hypers: "Psi" must be square')
      if (sum( !( (Hypers$T$Psi) == t(Hypers$T$Psi) ) ) > 0) stop('Hypers: "Psi" must be symmetric')
      if ((det(Hypers$T$Psi) - 0) < 0.00001) stop('Hypers: "Psi" is close to singular')
    }

    ###If phi hyperparameters are provided
    if ("Phi" %in% names(Hypers)) {
      if (!is.list(Hypers$Phi)) stop('Hypers: "Phi" must be a list')
      if (!"APhi" %in% names(Hypers$Phi)) stop('Hypers: "APhi" value missing')
      if (!is.scalar(Hypers$Phi$APhi)) stop('Hypers: "APhi" must be a scalar')
      if (is.na(Hypers$Phi$APhi)) stop('Hypers: "APhi" cannot be NA')
      if (!is.finite(Hypers$Phi$APhi)) stop('Hypers: "APhi" cannot be infinite')
      if (!"BPhi" %in% names(Hypers$Phi)) stop('Hypers: "BPhi" value missing')
      if (!is.scalar(Hypers$Phi$BPhi)) stop('Hypers: "BPhi" must be a scalar')
      if (is.na(Hypers$Phi$BPhi)) stop('Hypers: "BPhi" cannot be NA')
      if (!is.finite(Hypers$Phi$BPhi)) stop('Hypers: "BPhi" cannot be infinite')
      if (TemporalStructure == "exponential") {
        if (Hypers$Phi$APhi <= 0) stop('Hypers: For "exponential" correlation "APhi" must be strictly positive')
        if (Hypers$Phi$BPhi <= 0) stop('Hypers: For "exponential" correlation "BPhi" must be strictly positive')
        if (Hypers$Phi$BPhi < Hypers$Phi$APhi) stop('Hypers: "BPhi" must be greater than "APhi"')
      }
      if (TemporalStructure == "ar1") {
        if (Hypers$Phi$APhi < 0) stop('Hypers: For "ar1" correlation "APhi" must be in [0, 1]')
        if (Hypers$Phi$APhi > 1) stop('Hypers: For "ar1" correlation "APhi" must be in [0, 1]')
        if (Hypers$Phi$BPhi < 0) stop('Hypers: For "ar1" correlation "BPhi" must be in [0, 1]')
        if (Hypers$Phi$BPhi > 1) stop('Hypers: For "ar1" correlation "BPhi" must be in [0, 1]')
        if (Hypers$Phi$BPhi < Hypers$Phi$APhi) stop('Hypers: "BPhi" must be greater than "APhi"')
      }
    }

  ###End Hyperparameters
  }

  ###Starting Values
  if (!is.null(Starting)) {
    if (!is.list(Starting)) stop('Starting must be a list')
    if (!all(names(Starting) %in% c("Delta", "T", "Phi"))) stop('Starting: Can only contain objects with names "Delta", "T" and "Phi"')

    ###If delta starting values is provided
    if ("Delta" %in% names(Starting)) {
      if (!is.numeric(Starting$Delta)) stop('Starting: "Delta" must be a vector')
      if (length(Starting$Delta) != 3) stop('Starting: "Delta" must be length 3')
      if (!all(!is.na(Starting$Delta))) stop('Starting: "Delta" cannot have missing values')
      if (!all(is.finite(Starting$Delta))) stop('Starting: "Delta" cannot have infinite values')
    }

    ###If T starting values is provided
    if ("T" %in% names(Starting)) {
      if (!is.matrix(Starting$T)) stop('Starting: "T" must be a matrix')
      if (!dim(Starting$T)[1] == 3) stop('Starting: "T" must be 3 dimensional')
      if (!dim(Starting$T)[2] == 3) stop('Starting: "T" must be square')
      if (!all(!is.na(Starting$T))) stop('Starting: "T" cannot have missing values')
      if (!all(is.finite(Starting$T))) stop('Starting: "T" cannot have infinite values')
      if (sum( !( (Starting$T) == t(Starting$T) ) ) > 0) stop('Starting: "T" must be symmetric')
      if ((det(Starting$T) - 0) < 0.00001) stop('Starting: "T" is close to singular')
    }

    ###If phi starting values is provided
    if ("Phi" %in% names(Starting)) {
      if (!is.scalar(Starting$Phi)) stop('Starting: "Phi" must be a scalar')
      if (is.na(Starting$Phi)) stop('Starting: "Phi" cannot be NA')
      if (!is.finite(Starting$Phi)) stop('Starting: "Phi" cannot be infinite')
      # I make sure that Phi is in [APhi, BPhi] in CreateHyPara();
    }

  ###End Starting Values
  }

  ###Tuning Values
  if (!is.null(Tuning)) {
    if (!is.list(Tuning)) stop('Tuning must be a list')
    if (!all(names(Tuning) %in% c("Theta2", "Theta3", "Phi"))) stop('Tuning: Can only contain objects with names "Theta2", "Theta3" and "Phi"')

    ###If theta2 tuning values are provided
    if ("Theta2" %in% names(Tuning)) {
      if (!is.numeric(Tuning$Theta2)) stop('Tuning: "Theta2" must be a vector')
      if (length(Tuning$Theta2) != Nu) stop(paste0('Tuning: "Theta2" must have length ', Nu))
      if (!all(!is.na(Tuning$Theta2))) stop('Tuning: "Theta2" cannot have missing values')
      if (!all(is.finite(Tuning$Theta2))) stop('Tuning: "Theta2" cannot have infinite values')
      if (any(Tuning$Theta2 < 0)) stop('Tuning: "Theta2" must have non-negative components')
    }

    ###If theta3 tuning values are provided
    if ("Theta3" %in% names(Tuning)) {
      if (!"Theta3" %in% names(Tuning)) stop('Tuning: "Theta3" value missing')
      if (!is.numeric(Tuning$Theta3)) stop('Tuning: "Theta3" must be a vector')
      if (length(Tuning$Theta3) != Nu) stop(paste0('Tuning: "Theta3" must have length ', Nu))
      if (!all(!is.na(Tuning$Theta3))) stop('Tuning: "Theta3" cannot have missing values')
      if (!all(is.finite(Tuning$Theta3))) stop('Tuning: "Theta3" cannot have infinite values')
      if (any(Tuning$Theta3 < 0)) stop('Tuning: "Theta3" must have non-negative components')
    }

    ###If phi tuning value is provided
    if ("Phi" %in% names(Tuning)) {
      if (!is.scalar(Tuning$Phi)) stop('Tuning: "Phi" must be a scalar')
      if (is.na(Tuning$Phi)) stop('Tuning: "Phi" cannot be NA')
      if (!is.finite(Tuning$Phi)) stop('Tuning: "Phi" cannot be infinite')
      if (Tuning$Phi < 0) stop('Tuning: "Phi" must be non-negative')
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
