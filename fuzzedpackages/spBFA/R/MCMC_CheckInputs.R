CheckInputs <- function(formula, data, dist, time, K, L, trials, starting, hypers, tuning, mcmc, family, temporal.structure, spatial.structure, gamma.shrinkage, include.space, clustering) {
  
  ###Data dimensions
  N <- dim(data)[1]
  M <- dim(dist)[1]
  Nu <- length(time)
  O <- N / (M * Nu)
  
  ###Model specifications
  if (!is.logical(gamma.shrinkage)) stop('"gamma.shrinkage" must be a logical')
  if (!is.logical(include.space)) stop('"include.space" must be a logical')
  if (!is.logical(clustering)) stop('"clustering" must be a logical')
  if (!clustering & gamma.shrinkage) stop('gamma.shrinkage prior cannot be used without clustering')

  ###Family
  if ((length(family) != O) & (length(family) != 1)) stop(paste0('family: must have 1 or O = ', O, ' entries'))
  bool <- logical(length = length(family))
  for (i in 1:length(family)) bool[i] <- !(family[i] %in% c("normal", "probit", "tobit", "binomial"))
  if (any(bool)) stop('family: All entries must be one of "normal", "probit", "tobit", or "binomial"')
  
  ###Temporal correlation
  if (!temporal.structure %in% c("ar1", "exponential")) stop('temporal.structure: must be one of "ar1" or "exponential"')

  ###Spatial correlation
  if (!spatial.structure %in% c("discrete", "continuous")) stop('spatial.structure: must be one of "discrete" or "continuous"')
  
  ###Data and Formula
  if (!(class(formula) == "formula")) stop('"formula" must be of class formula')
  if (!is.data.frame(data)) stop('"data" must be of class data.frame')
  formula.test <- model.frame(formula, data) # check that the variables in formula are in data (will return error)
  
  ###Data checks for trials
  if (is.null(trials) & any(family %in% "binomial")) stop('trials must be specified for family "binomial"')
  if (!is.null(trials)) {
    Trials <- data[, trials] # checks if trials is in data
    if (!any(family %in% "binomial")) stop('trials can not be specified without family = "binomial"')
    if (length(Trials) != N) stop(paste0('trials must have exactly ', N, ' values')) # I check dimensions later
    if (any(is.na(Trials))) stop("trials may have no missing values")
    if (any(!is.finite(Trials))) stop("trials must have strictly finite entries")
    if (!isTRUE(all(Trials == floor(Trials)))) stop("trials must have integers only")
    if (any(Trials < 1)) stop("trials must contain positive integers only")
  }
  
  ###Data checks for Y
  familyInd <- numeric(length = O)
  if (length(family) == O) {
    for (o in 1:O) {
      if (family[o] == "normal") familyInd[o] <- "normal"
      if (family[o] == "probit") familyInd[o] <- "probit"
      if (family[o] == "tobit") familyInd[o] <- "tobit"
      if (family[o] == "binomial") familyInd[o] <- "binomial"
    }
  }
  if (length(family) == 1) {
    if (family == "normal") familyInd <- rep("normal", O)
    if (family == "probit") familyInd <- rep("probit", O)
    if (family == "tobit") familyInd <- rep("tobit", O)
    if (family == "binomial") familyInd <- rep("binomial", O)
  }
  outcome <- all.vars(formula)[1] # outcome variable
  Y <- formula.test[, outcome]
  YArray <- array(Y, dim = c(M, O, Nu))
  if (length(Y) != N) stop(paste0('Y must have exactly ', N, 'values'))
  if (any(is.na(Y))) stop("Y may have no missing values")
  if (any(!is.finite(Y))) stop("Y must have strictly finite entries")
  if (any(familyInd == "probit")) {
    for (o in 1:O) {
      if (familyInd[o] == "probit") {
        if ((sum(YArray[ , o, ] == 1) + sum(YArray[ , o, ] == 0)) != (Nu * M)) stop('Y: for "probit" observed data must be binary')
      }
    }
  } 
  if (any(familyInd == "tobit")) {
    for (o in 1:O) {
      if (familyInd[o] == "tobit") {
        if (any(YArray[ , o, ] < 0)) stop('Y: for "tobit" observed data must be non-negative')
      }
    }
  } 
  if (any(familyInd == "binomial")) {
    TrialsArray <- array(Trials, dim = c(M, O, Nu))
    count <- 1
    for (o in 1:O) {
      if (familyInd[o] == "binomial") {
        if (any(YArray[ , o, ] < 0)) stop('Y: for "binomial" observed data must be non-negative')
        if (!isTRUE(all(YArray[, o, ] == floor(YArray[ , o, ])))) stop('Y: for "binomial" observed data must be non-negative integers')
        if (any(YArray[, o, ] > TrialsArray[, count, ])) stop('Y: for "binomial" observed data must be less than the corresponding number of trials')
      }
    }
  } 
  
  ###Family indicator
  FamilyInd <- numeric(length = O)
  if (length(family) == O) {
    for (o in 1:O) {
      if (family[o] == "normal") FamilyInd[o] <- 0
      if (family[o] == "probit") FamilyInd[o] <- 1
      if (family[o] == "tobit") FamilyInd[o] <- 2
      if (family[o] == "binomial") FamilyInd[o] <- 3
    }
  }
  if (length(family) == 1) {
    if (family == "normal") FamilyInd <- rep(0, O)
    if (family == "probit") FamilyInd <- rep(1, O)
    if (family == "tobit") FamilyInd <- rep(2, O)
    if (family == "binomial") FamilyInd <- rep(3, O)
  }
  C <- sum(FamilyInd == 3)
  
  ###Data checks for covariates
  covariates <- all.vars(formula)[-1] # outcome variable
  X <- as.matrix(formula.test[, covariates])
  if (any(is.na(X))) stop('Covariates in formula cannot contain missing values')
  if (any(is.infinite(X))) stop('Covariates in formula cannot contain infinite values')
  
  ###Data checks for dist
  if (!is.matrix(dist)) stop('dist must be a matrix')
  if (!dim(dist)[1] == M) stop(paste0('dist must be a ', M ,' x ', M, ' dimensional matrix'))
  if (!dim(dist)[2] == M) stop('dist must be square')
  if (sum(!((dist) == t(dist))) > 0) stop('dist must be symmetric')
  if (any(diag(dist) != 0)) stop('dist must have only zeros on the diagonal')
  if (!all(!is.na(dist))) stop('dist cannot have missing values')
  if (!all(is.finite(dist))) stop('dist cannot have infinite values')
  if (spatial.structure == "discrete") if (length(table(dist)) > 2) stop('dist must only contain binaries (i.e. 0\'s or 1\'s) for "discrete" space')
  
  ###Data checks for time
  if (!is.numeric(time)) stop('time must be a vector')
  if (length(time) != Nu) stop(paste0('time must have length ', Nu))
  if (any(is.na(time))) stop("time may have no missing values")
  if (any(!is.finite(time))) stop("time must have strictly finite entries")
  if (is.unsorted(time)) stop('time vector is not in increasing order')
  if (!all(time >= 0)) stop('time vector has at least one negative point')

  ###Data checks for K
  if (missing(K)) stop("K: missing")
  if (!is.scalar(K)) stop('K must be a scalar')
  if (is.na(K)) stop('K cannot be NA')
  if (!is.finite(K)) stop('K cannot be infinite')
  if (!is.wholenumber(K) | K <= 0) stop('K must be a strictly positive integer')

  ###Data checks for L
  if (!is.scalar(L) & !is.infinite(L)) stop('L must be a scalar or Inf')
  if (is.na(L)) stop('L cannot be NA')
  if (is.finite(L)) if (!is.wholenumber(L) | L <= 0) stop('L as a scalar must be a strictly positive integer')
  
  ###hypers
  if (!is.null(hypers)) {
    if (!is.list(hypers)) stop('hypers must be a list')
    if (!all(names(hypers) %in% c("Sigma2", "Kappa", "Rho", "Delta", "Psi", "Upsilon"))) stop('hypers: Can only contain lists with names "Sigma2", "Kappa", "Rho", "Delta", "Psi", or "Upsilon"')

    ###If Sigma2 hyperparameters are provided
    if ("Sigma2" %in% names(hypers)) {
      if (any(family %in% c("normal", "probit", "tobit"))) {
        if (!is.list(hypers$Sigma2)) stop('hypers: "Sigma2" must be a list')
        if (!"A" %in% names(hypers$Sigma2)) stop('hypers: "A" value missing')
        if (!is.scalar(hypers$Sigma2$A)) stop('hypers: "A" must be a scalar')
        if (is.na(hypers$Sigma2$A)) stop('hypers: "A" cannot be NA')
        if (!is.finite(hypers$Sigma2$A)) stop('hypers: "A" cannot be infinite')
        if (hypers$Sigma2$A <= 0) stop('hypers: "A" must be strictly positive')
        if (!"B" %in% names(hypers$Sigma2)) stop('hypers: "B" value missing')
        if (!is.scalar(hypers$Sigma2$B)) stop('hypers: "B" must be a scalar')
        if (is.na(hypers$Sigma2$B)) stop('hypers: "B" cannot be NA')
        if (!is.finite(hypers$Sigma2$B)) stop('hypers: "B" cannot be infinite')
        if (hypers$Sigma2$B <= 0) stop('hypers: "B" must be strictly positive')
      } else stop('hypers: "Sigma2" cannot be included for "binomial" likelihood')
    }
    
    ###If Kappa hyperparameters are provided
    if ("Kappa" %in% names(hypers)) {
      if (!is.list(hypers$Kappa)) stop('hypers: "Kappa" must be a list')
      if (!"SmallUpsilon" %in% names(hypers$Kappa)) stop('hypers: "Kappa" value missing for SmallUpsilon')
      if (!is.scalar(hypers$Kappa$SmallUpsilon)) stop('hypers: "SmallUpsilon" must be a scalar')
      if (is.na(hypers$Kappa$SmallUpsilon)) stop('hypers: "SmallUpsilon" cannot be NA')
      if (!is.finite(hypers$Kappa$SmallUpsilon)) stop('hypers: "SmallUpsilon" cannot be infinite')
      if (hypers$Kappa$SmallUpsilon <= 0) stop('hypers: "SmallUpsilon" must be strictly positive')
      if (!"BigTheta" %in% names(hypers$Kappa)) stop('hypers: "BigTheta" value missing')
      if (O > 1) {
        if (hypers$Kappa$SmallUpsilon < O) stop('hypers: "SmallUpsilon" must be greater than or equal to O')
        if (!is.matrix(hypers$Kappa$BigTheta)) stop('hypers: "BigTheta" must be a matrix')
        if (!dim(hypers$Kappa$BigTheta)[1] == O) stop('hypers: "BigTheta" must be O dimensional')
        if (!all(!is.na(hypers$Kappa$BigTheta))) stop('hypers: "BigTheta" cannot have missing values')
        if (!all(is.finite(hypers$Kappa$BigTheta))) stop('hypers: "BigTheta" cannot have infinite values')
        if (!dim(hypers$Kappa$BigTheta)[2] == O) stop('hypers: "BigTheta" must be square')
        if (sum( !( (hypers$Kappa$BigTheta) == t(hypers$Kappa$BigTheta) ) ) > 0) stop('hypers: "BigTheta" must be symmetric')
        if ((det(hypers$Kappa$BigTheta) - 0) < 0.00001) stop('hypers: "BigTheta" is close to singular')
      }
      if (O == 1) {
        if (!is.scalar(hypers$Kappa$BigTheta)) stop('hypers: "BigTheta" must be a scalar when O is 1')
        if (is.na(hypers$Kappa$BigTheta)) stop('hypers: "BigTheta" cannot be NA')
        if (!is.finite(hypers$Kappa$BigTheta)) stop('hypers: "BigTheta" cannot be infinite')
        if (hypers$Kappa$BigTheta <= 0) stop('hypers: "BigTheta" must be strictly positive')
      }
    }

    ###If Rho hyperparameters are provided
    if ("Rho" %in% names(hypers)) {
      if (spatial.structure == "discrete") if (!is.null(hypers$Rho)) stop('hypers: When spatial.structure = "discrete", "Rho" must be missing or NULL')
      if (spatial.structure == "continuous") {
        if (!is.list(hypers$Rho)) stop('hypers: "Rho" must be a list')
        if (!"ARho" %in% names(hypers$Rho)) stop('hypers: "ARho" value missing')
        if (!is.scalar(hypers$Rho$ARho)) stop('hypers: "ARho" must be a scalar')
        if (is.na(hypers$Rho$ARho)) stop('hypers: "ARho" cannot be NA')
        if (!is.finite(hypers$Rho$ARho)) stop('hypers: "ARho" cannot be infinite')
        if (hypers$Rho$ARho <= 0) stop('hypers: "ARho" must be strictly positive')
        if (!"BRho" %in% names(hypers$Rho)) stop('hypers: "BRho" value missing')
        if (!is.scalar(hypers$Rho$BRho)) stop('hypers: "BRho" must be a scalar')
        if (is.na(hypers$Rho$BRho)) stop('hypers: "BRho" cannot be NA')
        if (!is.finite(hypers$Rho$BRho)) stop('hypers: "BRho" cannot be infinite')
        if (hypers$Rho$BRho <= 0) stop('hypers: "BRho" must be strictly positive')
        if (hypers$Rho$BRho < hypers$Rho$ARho) stop('hypers: "BRho" must be greater than "ARho"')
      }
    }
    
    ###If Delta hyperparameters are provided
    if ("Delta" %in% names(hypers)) {
      if (!is.list(hypers$Delta)) stop('hypers: "Delta" must be a list')
      if (!"A1" %in% names(hypers$Delta)) stop('hypers: "A1" value missing')
      if (!is.scalar(hypers$Delta$A1)) stop('hypers: "A1" must be a scalar')
      if (is.na(hypers$Delta$A1)) stop('hypers: "A1" cannot be NA')
      if (!is.finite(hypers$Delta$A1)) stop('hypers: "A1" cannot be infinite')
      if (hypers$Delta$A1 <= 0) stop('hypers: "A1" must be strictly positive')
      if (!"A2" %in% names(hypers$Delta)) stop('hypers: "A2" value missing')
      if (!is.scalar(hypers$Delta$A2)) stop('hypers: "A2" must be a scalar')
      if (is.na(hypers$Delta$A2)) stop('hypers: "A2" cannot be NA')
      if (!is.finite(hypers$Delta$A2)) stop('hypers: "A2" cannot be infinite')
      if (hypers$Delta$A2 <= 0) stop('hypers: "A2" must be strictly positive')
    }
    
    ###If Upsilon hyperparameters are provided
    if ("Upsilon" %in% names(hypers)) {
      if (!is.list(hypers$Upsilon)) stop('hypers: "Upsilon" must be a list')
      if (!"Zeta" %in% names(hypers$Upsilon)) stop('hypers: "Zeta" value missing')
      if (!is.scalar(hypers$Upsilon$Zeta)) stop('hypers: "Zeta" must be a scalar')
      if (is.na(hypers$Upsilon$Zeta)) stop('hypers: "Zeta" cannot be NA')
      if (!is.finite(hypers$Upsilon$Zeta)) stop('hypers: "Zeta" cannot be infinite')
      if (!"Omega" %in% names(hypers$Upsilon)) stop('hypers: "Omega" value missing')
      if (K > 1) {
        if (hypers$Upsilon$Zeta < K) stop('hypers: "Zeta" must be greater than or equal to K')
        if (!is.matrix(hypers$Upsilon$Omega)) stop('hypers: "Omega" must be a matrix')
        if (!dim(hypers$Upsilon$Omega)[1] == K) stop('hypers: "Omega" must be K dimensional')
        if (!all(!is.na(hypers$Upsilon$Omega))) stop('hypers: "Omega" cannot have missing values')
        if (!all(is.finite(hypers$Upsilon$Omega))) stop('hypers: "Omega" cannot have infinite values')
        if (!dim(hypers$Upsilon$Omega)[2] == K) stop('hypers: "Omega" must be square')
        if (sum( !( (hypers$Upsilon$Omega) == t(hypers$Upsilon$Omega) ) ) > 0) stop('hypers: "Omega" must be symmetric')
        if ((det(hypers$Upsilon$Omega) - 0) < 0.00001) stop('hypers: "Omega" is close to singular')
      }
      if (K == 1) {
        if (!is.scalar(hypers$Upsilon$Omega)) stop('hypers: "Omega" must be a scalar when K = 1')
        if (is.na(hypers$Upsilon$Omega)) stop('hypers: "Omega" cannot be NA')
        if (!is.finite(hypers$Upsilon$Omega)) stop('hypers: "Omega" cannot be infinite')
      }
    }

    ###If Psi hyperparameters are provided
    if ("Psi" %in% names(hypers)) {
      if (!is.list(hypers$Psi)) stop('hypers: "Psi" must be a list')
      if (temporal.structure == "exponential") {
        if (!"APsi" %in% names(hypers$Psi)) stop('hypers: "APsi" value missing')
        if (!is.scalar(hypers$Psi$APsi)) stop('hypers: "APsi" must be a scalar')
        if (is.na(hypers$Psi$APsi)) stop('hypers: "APsi" cannot be NA')
        if (!is.finite(hypers$Psi$APsi)) stop('hypers: "APsi" cannot be infinite')
        if (!"BPsi" %in% names(hypers$Psi)) stop('hypers: "BPsi" value missing')
        if (!is.scalar(hypers$Psi$BPsi)) stop('hypers: "BPsi" must be a scalar')
        if (is.na(hypers$Psi$BPsi)) stop('hypers: "BPsi" cannot be NA')
        if (!is.finite(hypers$Psi$BPsi)) stop('hypers: "BPsi" cannot be infinite')
        if (hypers$Psi$APsi < 0) stop('hypers: "APsi" must be non-negative')
        if (hypers$Psi$BPsi <= 0) stop('hypers: "BPsi" must be strictly positive')
        if (hypers$Psi$BPsi < hypers$Psi$APsi) stop('hypers: "BPsi" must be greater than "APsi"')
      }
      if (temporal.structure == "ar1") {
        if (!"APsi" %in% names(hypers$Psi)) stop('hypers: "APsi" value missing')
        if (!is.scalar(hypers$Psi$APsi)) stop('hypers: "APsi" must be a scalar')
        if (is.na(hypers$Psi$APsi)) stop('hypers: "APsi" cannot be NA')
        if (!is.finite(hypers$Psi$APsi)) stop('hypers: "APsi" cannot be infinite')
        if (!"BPsi" %in% names(hypers$Psi)) stop('hypers: "BPsi" value missing')
        if (!is.scalar(hypers$Psi$BPsi)) stop('hypers: "BPsi" must be a scalar')
        if (is.na(hypers$Psi$BPsi)) stop('hypers: "BPsi" cannot be NA')
        if (!is.finite(hypers$Psi$BPsi)) stop('hypers: "BPsi" cannot be infinite')
        if ((hypers$Psi$APsi < -1) | (hypers$Psi$APsi > 1)) stop('hypers: "APsi" must be in the range (-1, 1)')
        if ((hypers$Psi$BPsi < -1) | (hypers$Psi$BPsi > 1)) stop('hypers: "BPsi" must be in the range (-1, 1)')
        if (hypers$Psi$BPsi < hypers$Psi$APsi) stop('hypers: "BPsi" must be greater than "APsi"')
        if (!"Beta" %in% names(hypers$Psi)) stop('hypers: "Beta" value missing')
        if (!is.scalar(hypers$Psi$Beta)) stop('hypers: "Beta" must be a scalar')
        if (is.na(hypers$Psi$Beta)) stop('hypers: "Beta" cannot be NA')
        if (!is.finite(hypers$Psi$Beta)) stop('hypers: "Beta" cannot be infinite')
        if (!"Gamma" %in% names(hypers$Psi)) stop('hypers: "Gamma" value missing')
        if (!is.scalar(hypers$Psi$Gamma)) stop('hypers: "Gamma" must be a scalar')
        if (is.na(hypers$Psi$Gamma)) stop('hypers: "Gamma" cannot be NA')
        if (!is.finite(hypers$Psi$Gamma)) stop('hypers: "Gamma" cannot be infinite')
        if (hypers$Psi$Beta <= 0) stop('hypers: "Beta" must be strictly positive')
        if (hypers$Psi$Gamma <= 0) stop('hypers: "Gamma" must be strictly positive')
      }
    }

  ###End Hyperparameters
  }

  ###starting Values
  if (!is.null(starting)) {
    if (!is.list(starting)) stop('starting must be a list')
    if (!all(names(starting) %in% c("Sigma2", "Kappa", "Rho", "Delta", "Psi", "Upsilon"))) stop('starting: Can only contain objects with names "Sigma2", "Kappa", "Rho", "Delta", "Psi", and "Upsilon"')

    ###If Delta starting values is provided
    if ("Delta" %in% names(starting)) {
      if (is.vector(starting$Delta)) {
        if (!is.numeric(starting$Delta)) stop('starting: "Delta" must be a vector')
        if (length(starting$Delta) != K) stop('starting: "Delta" must be length K')
        if (!all(!is.na(starting$Delta))) stop('starting: "Delta" cannot have missing values')
        if (!all(is.finite(starting$Delta))) stop('starting: "Delta" cannot have infinite values')
        if (any(starting$Delta < 0)) stop('starting: "Delta" cannot have non-negative values')
      }
    }
    if (is.scalar(starting$Delta)) {
      if (is.na(starting$Delta)) stop('starting: "Delta" cannot be NA')
      if (!is.finite(starting$Delta)) stop('starting: "Delta" cannot be infinite')
      if (starting$Delta < 0) stop('starting: "Delta" must be non-negative')
    }
    if ((!is.vector(starting$Delta)) & (!is.scalar(starting$Delta))) stop('starting: "Delta" must be a scalar or a vector')
    
    ###If Sigma2 starting values is provided
    if ("Sigma2" %in% names(starting)) {
      if (any(family %in% c("normal", "probit", "tobit"))) {
        if (is.matrix(starting$Sigma2)) {
          if (!dim(starting$Sigma2)[1] == M) stop(paste0('starting: "Sigma2" must have', M, 'rows'))
          if (!dim(starting$Sigma2)[2] == (O - C)) stop(paste0('starting: "Sigma2" must have', O - C, 'columns'))
          if (!all(!is.na(starting$Sigma2))) stop('starting: "Sigma2" cannot have missing values')
          if (!all(is.finite(starting$Sigma2))) stop('starting: "Sigma2" cannot have infinite values')
          if (starting$Sigma2 <= 0) stop('starting: "Sigma2" must be strictly positive')
        }
        if (is.scalar(starting$Sigma2)) {
            if (is.na(starting$Sigma2)) stop('starting: "Sigma2" cannot be NA')
            if (!is.finite(starting$Sigma2)) stop('starting: "Sigma2" cannot be infinite')
            if (starting$Sigma2 < 0) stop('starting: "Sigma2" must be non-negative')
        }
        if ((!is.matrix(starting$Sigma2)) & (!is.scalar(starting$Sigma2))) stop('starting: "Sigma2" must be a scalar or a matrix')
      } else stop('starting: "Sigma2" does not get included for "binomial" likelihood')
    }
    
    ###If Kappa starting values is provided
    if ("Kappa" %in% names(starting)) {
      if (O > 1) {
        if (!is.matrix(starting$Kappa)) stop('starting: "Kappa" must be a matrix')
        if (!dim(starting$Kappa)[1] == O) stop('starting: "Kappa" must be O dimensional')
        if (!dim(starting$Kappa)[2] == O) stop('starting: "Kappa" must be square')
        if (!all(!is.na(starting$Kappa))) stop('starting: "Kappa" cannot have missing values')
        if (!all(is.finite(starting$Kappa))) stop('starting: "Kappa" cannot have infinite values')
        if (sum( !( (starting$Kappa) == t(starting$Kappa) ) ) > 0) stop('starting: "Kappa" must be symmetric')
        if ((det(starting$Kappa) - 0) < 0.0000000001) stop('starting: "Kappa" is close to singular')
      }
      if (O == 1) {
        if (!is.scalar(starting$Kappa)) stop('starting: "Kappa" must be a scalar when O = 1')
        if (is.na(starting$Kappa)) stop('starting: "Kappa" cannot be NA')
        if (!is.finite(starting$Kappa)) stop('starting: "Kappa" cannot be infinite')
        if (starting$Kappa <= 0) stop('starting: "Kappa" must be positive')
      }
    }

    ###If Rho starting values is provided
    if ("Rho" %in% names(starting)) {
      if (!is.scalar(starting$Rho)) stop('starting: "Rho" must be a scalar')
      if (is.na(starting$Rho)) stop('starting: "Rho" cannot be NA')
      if (!is.finite(starting$Rho)) stop('starting: "Rho" cannot be infinite')
      if (spatial.structure == "discrete") {
        if ((starting$Rho <= 0) | (starting$Rho >= 1)) stop('starting: "Rho" must be in (0, 1) for discrete spatial process')
      }
      if (spatial.structure == "continuous") {
        # I make sure that Rho is in (ARho, BRho) in CreatePara();
      }
    }
    
    ###If Upsilon starting values is provided
    if ("Upsilon" %in% names(starting)) {
      if (K > 1) {
        if (!is.matrix(starting$Upsilon)) stop('starting: "Upsilon" must be a matrix')
        if (!dim(starting$Upsilon)[1] == K) stop('starting: "Upsilon" must be K dimensional')
        if (!dim(starting$Upsilon)[2] == K) stop('starting: "Upsilon" must be square')
        if (!all(!is.na(starting$Upsilon))) stop('starting: "Upsilon" cannot have missing values')
        if (!all(is.finite(starting$Upsilon))) stop('starting: "Upsilon" cannot have infinite values')
        if (sum( !( (starting$Upsilon) == t(starting$Upsilon) ) ) > 0) stop('starting: "Upsilon" must be symmetric')
        if ((det(starting$Upsilon) - 0) < 0.0000000001) stop('starting: "Upsilon" is close to singular')
      }
      if (K == 1) {
        if (!is.scalar(starting$Upsilon)) stop('starting: "Upsilon" must be a scalar when K = 1')
        if (is.na(starting$Upsilon)) stop('starting: "Upsilon" cannot be NA')
        if (!is.finite(starting$Upsilon)) stop('starting: "Upsilon" cannot be infinite')
        if (starting$Upsilon <= 0) stop('starting: "Upsilon" must be positive')
      }
    }

    ###If Psi starting values is provided
    if ("Psi" %in% names(starting)) {
      if (!is.scalar(starting$Psi)) stop('starting: "Psi" must be a scalar')
      if (is.na(starting$Psi)) stop('starting: "Psi" cannot be NA')
      if (!is.finite(starting$Psi)) stop('starting: "Psi" cannot be infinite')
      # I make sure that Psi is in (APsi, BPsi) in CreatePara();
    }

  ###End starting Values
  }

  ###tuning Values
  if (!is.null(tuning)) {
    if (!is.list(tuning)) stop('tuning must be a list')
    if (!all(names(tuning) %in% c("Psi", "Rho"))) stop('tuning: Can only contain objects with names "Psi", "Rho"')

    ###If Psi tuning value is provided
    if ("Psi" %in% names(tuning)) {
      if (!is.scalar(tuning$Psi)) stop('tuning: "Psi" must be a scalar')
      if (is.na(tuning$Psi)) stop('tuning: "Psi" cannot be NA')
      if (!is.finite(tuning$Psi)) stop('tuning: "Psi" cannot be infinite')
      if (tuning$Psi < 0) stop('tuning: "Psi" must be non-negative')
    }

    ###If Rho tuning value is provided
    if ("Rho" %in% names(tuning)) {
      if (spatial.structure == "discrete") {
        if (!is.null(tuning$Rho)) stop('tuning: No tuning is needed for "Rho" when "spatial.structure" is discrete')
      }
      if (spatial.structure == "continuous") {
        if (!is.scalar(tuning$Rho)) stop('tuning: "Rho" must be a scalar')
        if (is.na(tuning$Rho)) stop('tuning: "Rho" cannot be NA')
        if (!is.finite(tuning$Rho)) stop('tuning: "Rho" cannot be infinite')
        if (tuning$Rho < 0) stop('tuning: "Rho" must be non-negative')
      }
    }

  ###End tuning Values
  }

  ###mcmc Values
  if (!is.null(mcmc)) {
    if (!is.list(mcmc)) stop('mcmc must be a list')
    if (!all(names(mcmc) %in% c("NBurn", "NSims", "NThin", "NPilot"))) stop('mcmc: Can only contain objects with names "NBurn", "NSims", "NThin" and "NPilot"')

    ###If NBurn is provided
    if ("NBurn" %in% names(mcmc)) {
      if (!is.scalar(mcmc$NBurn)) stop('mcmc: "NBurn" must be a scalar')
      if (is.na(mcmc$NBurn)) stop('mcmc: "NBurn" cannot be NA')
      if (!is.finite(mcmc$NBurn)) stop('mcmc: "NBurn" cannot be infinite')
      if (!is.wholenumber(mcmc$NBurn) | mcmc$NBurn < 0) stop('mcmc: "NBurn" must be a non-negative integer')
      if (mcmc$NBurn < 100) stop('mcmc: "NBurn" must be at least 100')
    }

    ###If NSims is provided
    if ("NSims" %in% names(mcmc)) {
      if (!is.scalar(mcmc$NSims)) stop('mcmc: "NSims" must be a scalar')
      if (is.na(mcmc$NSims)) stop('mcmc: "NSims" cannot be NA')
      if (!is.finite(mcmc$NSims)) stop('mcmc: "NSims" cannot be infinite')
      if (!is.wholenumber(mcmc$NSims) | mcmc$NSims <= 0) stop('mcmc: "NSims" must be a positive integer')
      if (mcmc$NSims < 100) stop('mcmc: "NSims" must be at least 100')
    }

    ###If NThin is provided
    if ("NThin" %in% names(mcmc)) {
      if (!is.scalar(mcmc$NThin)) stop('mcmc: "NThin" must be a scalar')
      if (is.na(mcmc$NThin)) stop('mcmc: "NThin" cannot be NA')
      if (!is.finite(mcmc$NThin)) stop('mcmc: "NThin" cannot be infinite')
      if (!is.wholenumber(mcmc$NThin) | mcmc$NThin <= 0) stop('mcmc: "NThin" must be a positive integer')
      # if (!is.wholenumber(mcmc$NSims / mcmc$NThin)) stop('mcmc: "NThin" must be a factor of "NSims"') enforced in Createmcmc();
    }

    ###If NPilot is provided
    if ("NPilot" %in% names(mcmc)) {
      if (!is.scalar(mcmc$NPilot)) stop('mcmc: "NPilot" must be a scalar')
      if (is.na(mcmc$NPilot)) stop('mcmc: "NPilot" cannot be NA')
      if (!is.finite(mcmc$NPilot)) stop('mcmc: "NPilot" cannot be infinite')
      if (!is.wholenumber(mcmc$NPilot) | mcmc$NPilot < 0) stop('mcmc: "NPilot" must be a positive integer')
      # if (!is.wholenumber(mcmc$NBurn / mcmc$NPilot)) stop('mcmc: "NPilot" must be a factor of "NBurn"') enforced in Createmcmc();
    }

  ###End mcmc Values
  }

}

###Helper Functions
is.scalar <- function(x) ((is.numeric(x)) & (length(x) == 1))
is.wholenumber <- function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol
