# Copyright 2015-2019 Venelin Mitov
#
# This file is part of POUMM.
#
# POUMM is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# POUMM is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with POUMM  If not, see <http://www.gnu.org/licenses/>.


#' @name specPOUMM
#' 
#' @title Specifying a POUMM fit
#' 
#' @description Specification and validation of POUMM/PMM settings.
#' 
#' @param z,tree a numeric vector and a phylo object on which the fit is to be done. 
#'    These arguments are used in order to guess meaningful values for the parLower,
#'    parUpper and parPriorMCMC arguments. See also, zMin,zMean,...,tMax below.
#' @param zMin,zMean,zMax,zVar,zSD,tMin,tMean,tMax summary statistics of the
#'   observed tip-values (z) and root-tip distances (t). Some of these values
#'    are used for constructing default parameter values and limits; These 
#'    arguments are given default values which will most likely be meaningless
#'    in your specific use-case. The default values will be overwritten with the
#'    corresponding statistics from the z and tree arguments if these were specified.
#'    If none of tree and z, nor these parameters are specified, then the
#'    arguments parLower, parUpper, parPriorMCMC must be specified explicitly.
#' @param parMapping An R-function that can handle, both, a numeric vector 
#'   or a numeric matrix as argument. This function should transform the input 
#'   vector or each row-vector (if the input is matrix) into a (row-)vector of 
#'   the POUMM parameters alpha, theta, sigma, sigmae, g0. For a vector input 
#'   the function should return a vector with named elements alpha, theta, 
#'   sigma, sigmae, g0. For a matrix input the function should return a matrix 
#'   with the same number of rows and columns alpha, theta, sigma, sigmae, g0. 
#'   Only finite non-negative values are allowed for alpha, sigma, and sigmae. 
#'   Returning Inf, -Inf, NA or NaN for any of these parameters will result in 
#'   an error during likelihood calculation. Only finite numerical values are 
#'   allowed for theta. The parameter
#'   g0 is treated in a special way and can assume either a finite numerical 
#'   value or one of NA or NaN. If g0 = finite value, this value is used
#'   together
#'   with the corresponding values of alpha, theta, sigma, and sigmae for 
#'   likelihood calcuation. If g0 = NA (meaing value Not Avaiable), the value of
#'   g0 is calculated analytically during likelihood calculation in order to 
#'   maximise one of the following: \cr
#'   \enumerate{
#'    \item if a normal prior for g0 was specified (see g0Prior), 
#'       \eqn{pdf(z | \alpha, \theta, \sigma, \sigma_e, g0, tree) x prior(g0)}. 
#'    \item otherwise, \eqn{pdf(z | \alpha, \theta, \sigma, \sigma_e, g0, tree)}. 
#'   }
#'   If g0 = NaN (meaning Not a Number), then the likelihood is marginalized w.r.t. 
#'    the g0's prior distribution (see g0Prior), i.e. the likelihood returned is:
#'    \eqn{pdf(z | \alpha, \theta, \sigma, \sigma_e, tree) = Integral(pdf(z|\alpha,\theta,\sigma,\sigma_e,g0) x pdf(g0) d g0; g0 from -\infty to +\infty) }
#'    In this case (g0=NaN), if g0Prior is not specified, 
#'    it is assumed that g0Prior is the stationary OU normal distribution with 
#'    mean, theta, and variance, varOU(Inf, alpha, sigma). \cr
#'  Examples: \cr
#'   
#'  \preformatted{ 
#'  # Default for POUMM: identity for alpha, theta, sigma, sigmae, NA for g0.
#'  parMapping = function(par) {
#'    if(is.matrix(par)) {
#'      atsseg0 <- cbind(par[, 1:4, drop = FALSE], NA) 
#'      colnames(atsseg0) <- c("alpha", "theta", "sigma", "sigmae", "g0")
#'    } else {
#'      atsseg0 <- c(par[1:4], NA) 
#'      names(atsseg0) <- c("alpha", "theta", "sigma", "sigmae", "g0")
#'    }
#'    atsseg0
#'  }
#' }
#' @param parLower,parUpper two named numeric vectors of the same length
#'   indicating the boundaries of the search region for the ML-fit. Calling
#'   parMapping on parLower and parUpper should result in appropriate values of
#'   the POUMM parameters alpha, theta, sigma sigmae and g0. By default, the
#'   upper limit for alpha is set to 69.31 / tMean, which corresponds to a value
#'   of alpha so big that the time for half-way convergence towards theta from
#'   any initial trait value is 100 times shorter than the mean root-tip distance
#'   in the tree. Examples: \cr
#' \preformatted{
#' # Default for POUMM:
#' parLower = c(alpha = 0, theta = zMin - 2 * (zMax - zMin), sigma = 0, sigmae = 0)
#' parUpper = c(alpha = 69.31 / tMean, theta = zMax + 2 * (zMax - zMin), 
#'              sigma = sigmaOU(H2 = .99, alpha = 69.31 / tMean, sigmae = 2 * zSD,
#'                                     t = tMean), 
#'              sigmae = 2 * zSD)
#' }
#' 
#' @param g0Prior Either NULL or a list with named numeric or character 
#'    members "mean" and "var". Specifies a prior normal distribution for the
#'    parameter g0. If characters, the members "mean" and "var" are evaluated as
#'    R-expressions - useful if these are functions of some of other parameters.
#'    Note that if g0Prior is not NULL and g0 is not NaN (either a fixed number
#'    or NA), then the likelihood maximization takes into account the prior for 
#'    g0, that is, the optimization is done over the product 
#'    p(g0) x lik(data|g0, other parameters and tree). This can be helpful to 
#'    prevent extremely big or low estimates of g0. To avoid this behavior and
#'    always maximize the likelihood, use g0Prior = NULL. 
#' @param sigmaeFixed fixed value for the sigmae parameter (used in 
#' specifyPOUMM_ATS and specifyPOUMM_ATSG0).
#' @param parInitML A named vector (like parLower and parUpper) or a list of such
#'   vectors - starting points for optim.
#' @param control List of parameters passed on to optim in the ML-fit, default 
#'   list(factr=1e9), see ?optim.
#'
#' @param parInitMCMC a function(chainNo, fitML) returning an initial state of
#'   an MCMC as a vector. The argument fitML can be used to specify an initial 
#'   state, close to a previously found likelihood optimum. Example: \cr
#' \preformatted{ 
#'  # Default for POUMM:
#'  parInitMCMC = function(chainNo, fitML) {
#'    if(!is.null(fitML)) {
#'      parML <- fitML$par
#'    } else {
#'      parML <- NULL
#'    }
#'    
#'    init <- rbind(
#'      c(alpha = 0, theta = 0, sigma = 1, sigmae = 0),
#'      parML,
#'      c(alpha = 0, theta = 0, sigma = 1, sigmae = 1)
#'    )
#'    
#'    init[(chainNo - 1) \%\% nrow(init) + 1, ]
#'  }
#' }
#' 
#' @param parPriorMCMC A function of a numeric parameter-vector returning the 
#' log-prior for this parameter vector. Example: \cr
#' 
#' \preformatted{
#' # Default for POUMM:
#'  parPriorMCMC = function(par) {
#'    dexp(par[1], rate = tMean / 6.931, TRUE) + 
#'      dnorm(par[2], zMean, 10 * zSD, TRUE) +
#'      dexp(par[3],  rate = sqrt(tMean / (zVar * 0.6931)), TRUE) + 
#'      dexp(par[4], rate = 2 / zSD, TRUE)
#'  }
#' }
#' 
#' @param parScaleMCMC Numeric matrix indicating the initial jump-distribution 
#'   matrix for the MCMC fit. Default for POUMM is diag(4); 
#' @param nSamplesMCMC Integer indicating the length of each MCMC chain. 
#' Defaults to 1e5.
#' @param nAdaptMCMC Logical indicating whether adaptation of the MCMC jump 
#'   distribution should be done with respect to the target acceptance rate 
#'   (accRateMCMC) or integer indicating how many initial MCMC iterations should 
#'   be used for adaptation of the jump-distribution matrix (see details in 
#'   ?POUMM). Defaults to nSamplesMCMC meaning continuous adaptation throughout
#'   the MCMC. 
#' @param thinMCMC Integer indicating the thinning interval of the mcmc-chains. 
#'   Defaults to 100.
#' @param accRateMCMC numeric between 0 and 1 indicating the target 
#'   acceptance rate of the  adaptive Metropolis sampling (see details in ?POUMM). 
#'   Default 0.01.
#' @param gammaMCMC controls the speed of adaption. Should be in the interval (0.5,1]. A lower gamma leads to faster adaption. Default value is 0.50001.
#' @param nChainsMCMC integer indicating the number of chains to run. 
#'   Defaults to 3 chains, from which the first one is a sample from the prior 
#'   distribution (see samplePriorMCMC).
#' @param samplePriorMCMC Logical indicating if sampling from the prior 
#'   should be done for the first chain (see nChainsMCMC). This is useful to 
#'   compare mcmc's for an overlap between prior and posterior distributions. 
#'   Default is TRUE.
#' @param parallelMCMC Logical indicating whether the MCMC chains should be run 
#'   in parallel. Setting this option to TRUE results in using 
#'   \code{foreach::foreach() \%dopar\% { }} construct for the MCMC fit. In order for
#'   parallel execution to be done, you should create a computing cluster and
#'   register it as parallel back-end (see example in package vignette and the 
#'   web-page https://github.com/tobigithub/R-parallel/wiki/R-parallel-Setups).
#' @param validateSpec Logical indicating whether the passed parameters should 
#'   be validated. This parameter is used internally and should always be TRUE.
#' @return A named list to be passed as a spec argument to POUMM.
NULL

#' @describeIn specPOUMM Specify parameters for fitting a POUMM model. 
#'   Parameter vector is c(alpha, theta, sigma, sigmae). Default model settings. 
#' 
#' @export
specifyPOUMM <- function(
  z = NULL, tree = NULL,
  zMin = -10, zMean = 0, zMax = 10, zVar = 4, zSD = sqrt(zVar), 
  tMin = 0.1, tMean = 2, tMax = 10, 
  parMapping = NULL, 
  parLower = NULL, parUpper = NULL, 
  g0Prior = NULL,
  parInitML = NULL,
  control = NULL,
  
  parPriorMCMC = NULL, 
  parInitMCMC = NULL, 
  parScaleMCMC = NULL,
  nSamplesMCMC = 1e5, nAdaptMCMC = nSamplesMCMC, 
  thinMCMC = 100, 
  accRateMCMC = .01, gammaMCMC = 0.50001, nChainsMCMC = 3, 
  samplePriorMCMC = TRUE,
  parallelMCMC = FALSE,
  validateSpec=TRUE) {
  
  spec <- list(parMapping = parMapping, 
               parLower = parLower, parUpper = parUpper, 
               g0Prior = g0Prior,
               parInitML = parInitML,
               control = control,
               
               parPriorMCMC = parPriorMCMC, 
               parInitMCMC = parInitMCMC, 
               parScaleMCMC = parScaleMCMC,
               nSamplesMCMC = nSamplesMCMC, nAdaptMCMC = nAdaptMCMC, 
               thinMCMC = thinMCMC, accRateMCMC = accRateMCMC, 
               gammaMCMC = gammaMCMC, nChainsMCMC = nChainsMCMC, 
               samplePriorMCMC = samplePriorMCMC,
               parallelMCMC = parallelMCMC)
  
  if(validateZTree(z, tree)) {
    zMin <- min(z); zMean <- mean(z); zMax <- max(z); 
    zVar <- var(z); zSD <- sd(z);
    tipTimes <- nodeTimes(tree, tipsOnly = TRUE)
    tMin <- min(tipTimes); tMean <- mean(tipTimes); tMax <- max(tipTimes); 
  }
  
  specDefault <- list(
    parMapping = function(par) {
      if(is.matrix(par)) {
        atsseg0 <- cbind(par[, 1:4, drop = FALSE], NA) 
        colnames(atsseg0) <- c("alpha", "theta", "sigma", "sigmae", "g0")
      } else {
        atsseg0 <- c(par[1:4], NA) 
        names(atsseg0) <- c("alpha", "theta", "sigma", "sigmae", "g0")
      }
      atsseg0
    },
    
    parLower = c(alpha = 0, theta = zMin - 2 * (zMax - zMin), 
                 sigma = 0, sigmae = 0), 
    
    parUpper = c(alpha = 69.31 / tMean, theta = zMax + 2 * (zMax - zMin), 
                 sigma = sigmaOU(H2 = .99, alpha = 69.31 / tMean, sigmae = 2 * zSD, t = tMean), 
                 sigmae = 2 * zSD),
    
    g0Prior = NULL, 
    
    parInitML = NULL,
    
    control = list(factr = 1e8),
    
    parPriorMCMC = function(par) {
      dexp(par[1], rate = tMean / 6.931, log = TRUE) + 
        dnorm(par[2], zMean, 2 * zSD, log = TRUE) +
        dexp(par[3],  rate = sqrt(tMean / (zVar * 0.6931)), log = TRUE) + 
        dexp(par[4], rate = 2 / zSD, log = TRUE)  
      
    },
    
    parInitMCMC = function(chainNo, fitML = NULL) {
      if(!is.null(fitML)) {
        parML <- fitML$par
      } else {
        parML <- NULL
      }
      
      init <- rbind(
        c(alpha = 0, theta = zMean, sigma = 1, sigmae = 0),
        parML,
        c(alpha = 0, theta = zMean, sigma = 1, sigmae = 1)
      )
      
      init[(chainNo - 1) %% nrow(init) + 1, ]
    },
    
    parScaleMCMC = diag(4),
    
    nSamplesMCMC = nSamplesMCMC, nAdaptMCMC = nAdaptMCMC, 
    thinMCMC = thinMCMC, accRateMCMC = accRateMCMC, gammaMCMC = gammaMCMC, 
    
    nChainsMCMC = nChainsMCMC, samplePriorMCMC = samplePriorMCMC
  )
  
  if(is.list(spec)) {
    for(name in names(specDefault)) {
      if(is.null(spec[[name]]) & !is.null(specDefault[[name]])) {
        spec[[name]] <- specDefault[[name]]
      } 
    } 
  } else {
    spec <- specDefault
  }
  
  if(validateSpec) {
    validateSpecPOUMM(spec)
  }
  
  spec
}

#' @describeIn specPOUMM Fitting a POU model with fixed sigmae.
#'  Parameter vector is c(alpha, theta, sigma).
#' @export
specifyPOUMM_ATS <- function(
  z = NULL, tree = NULL,
  zMin = -10, zMean = 0, zMax = 10, zVar = 4, zSD = sqrt(zVar), 
  tMin = 0.1, tMean = 2, tMax = 10, 
  parMapping = NULL, 
  parLower = NULL, parUpper = NULL, 
  g0Prior = NULL,
  parInitML = NULL,
  control = NULL,
  parPriorMCMC = NULL, 
  parInitMCMC = NULL, 
  parScaleMCMC = NULL,
  nSamplesMCMC = 1e5, nAdaptMCMC = nSamplesMCMC, 
  thinMCMC = 100, 
  accRateMCMC = .01, gammaMCMC = 0.50001, nChainsMCMC = 3, 
  samplePriorMCMC = TRUE,
  parallelMCMC = FALSE, 
  sigmaeFixed = 0) {
  
  spec <- list(parMapping = parMapping, 
               parLower = parLower, parUpper = parUpper, 
               g0Prior = g0Prior,
               parInitML = parInitML,
               control = control,
               parPriorMCMC = parPriorMCMC, 
               parInitMCMC = parInitMCMC, 
               parScaleMCMC = parScaleMCMC,
               nSamplesMCMC = nSamplesMCMC, nAdaptMCMC = nAdaptMCMC, 
               thinMCMC = thinMCMC, accRateMCMC = accRateMCMC, 
               gammaMCMC = gammaMCMC, nChainsMCMC = nChainsMCMC, 
               samplePriorMCMC = samplePriorMCMC,
               parallelMCMC = parallelMCMC)
  
  if(validateZTree(z, tree)) {
    zMin <- min(z); zMean <- mean(z); zMax <- max(z); 
    zVar <- var(z); zSD <- sd(z);
    tipTimes <- nodeTimes(tree, tipsOnly = TRUE)
    tMin <- min(tipTimes); tMean <- mean(tipTimes); tMax <- max(tipTimes); 
  }
  
  specDefault <- specifyPOUMM(
    z, tree, zMin, zMean, zMax, zVar, zSD, tMin, tMean, tMax, 
    
    parMapping = function(par) {
      if(is.matrix(par)) {
        atsseg0 <- cbind(par[, 1:3, drop = FALSE], sigmaeFixed, NA) 
        colnames(atsseg0) <- c("alpha", "theta", "sigma", "sigmae", "g0")
      } else {
        atsseg0 <- c(par[1:3], sigmaeFixed, NA) 
        names(atsseg0) <- c("alpha", "theta", "sigma", "sigmae", "g0")
      }
      atsseg0
    },
    
    parLower = c(alpha = 0, theta = zMin - 2 * (zMax - zMin), 
                 sigma = 0), 
    
    parUpper = c(alpha = 69.31 / tMean, theta = zMax + 2 * (zMax - zMin), 
                 sigma = sigmaOU(H2 = .99, alpha = 69.31 / tMean, sigmae = 2 * zSD, t = tMean)),
    
    parPriorMCMC = function(par) {
      dexp(par[1], rate = tMean / 6.931, log = TRUE) +
        dnorm(par[2], zMean, 2 * zSD, log = TRUE) +
        dexp(par[3],  rate = sqrt(tMean / (zVar * 0.6931)), log = TRUE)
    },
    
    parInitMCMC = function(chainNo, fitML = NULL) {
      if(!is.null(fitML)) {
        parML <- fitML$par
      } else {
        parML <- NULL
      }
      
      init <- rbind(
        c(alpha = 0, theta = zMean, sigma = 1),
        parML,
        c(alpha = 0, theta = zMean, sigma = 1)
      )
      
      init[(chainNo - 1) %% nrow(init) + 1, ]
    },
    
    parScaleMCMC = diag(3),
    
    validateSpec = FALSE
  )
  
  if(is.list(spec)) {
    for(name in names(specDefault)) {
      if(is.null(spec[[name]]) & !is.null(specDefault[[name]])) {
        spec[[name]] <- specDefault[[name]]
      } 
    } 
  } else {
    spec <- specDefault
  }
  
  validateSpecPOUMM(spec)
  spec
}

#' @describeIn specPOUMM Fitting a POU model with fixed sigmae.
#'  Parameter vector is c(alpha, theta, sigma, g0).
#' @export
specifyPOUMM_ATSG0 <- function(
  z = NULL, tree = NULL,
  zMin = -10, zMean = 0, zMax = 10, zVar = 4, zSD = sqrt(zVar), 
  tMin = 0.1, tMean = 2, tMax = 10, 
  parMapping = NULL, 
  parLower = NULL, parUpper = NULL, 
  g0Prior = NULL,
  parInitML = NULL,
  control = NULL,
  parPriorMCMC = NULL, 
  parInitMCMC = NULL, 
  parScaleMCMC = NULL,
  nSamplesMCMC = 1e5, nAdaptMCMC = nSamplesMCMC, 
  thinMCMC = 100, 
  accRateMCMC = .01, gammaMCMC = 0.50001, nChainsMCMC = 3, 
  samplePriorMCMC = TRUE,
  parallelMCMC = FALSE, 
  sigmaeFixed = 0) {
  
  spec <- list(parMapping = parMapping, 
               parLower = parLower, parUpper = parUpper, 
               g0Prior = g0Prior,
               parInitML = parInitML,
               control = control,
               parPriorMCMC = parPriorMCMC, 
               parInitMCMC = parInitMCMC, 
               parScaleMCMC = parScaleMCMC,
               nSamplesMCMC = nSamplesMCMC, nAdaptMCMC = nAdaptMCMC, 
               thinMCMC = thinMCMC, accRateMCMC = accRateMCMC, 
               gammaMCMC = gammaMCMC, nChainsMCMC = nChainsMCMC, 
               samplePriorMCMC = samplePriorMCMC,
               parallelMCMC = parallelMCMC)
  
  if(validateZTree(z, tree)) {
    zMin <- min(z); zMean <- mean(z); zMax <- max(z); 
    zVar <- var(z); zSD <- sd(z);
    tipTimes <- nodeTimes(tree, tipsOnly = TRUE)
    tMin <- min(tipTimes); tMean <- mean(tipTimes); tMax <- max(tipTimes); 
  }
  
  specDefault <- specifyPOUMM(
    z, tree, zMin, zMean, zMax, zVar, zSD, tMin, tMean, tMax, 
    
    parMapping = function(par) {
      if(is.matrix(par)) {
        atsseg0 <- cbind(par[, 1:3, drop = FALSE], sigmaeFixed, par[, 4, drop = FALSE]) 
        colnames(atsseg0) <- c("alpha", "theta", "sigma", "sigmae", "g0")
      } else {
        atsseg0 <- c(par[1:3], sigmaeFixed, par[4]) 
        names(atsseg0) <- c("alpha", "theta", "sigma", "sigmae", "g0")
      }
      atsseg0
    },
    
    parLower = c(alpha = 0, theta = zMin - 2 * (zMax - zMin), 
                 sigma = 0, 
                 g0 = zMin - 2 * zSD), 
    
    parUpper = c(alpha = 69.31 / tMean, theta = zMax + 2 * (zMax - zMin), 
                 sigma = sigmaOU(H2 = .99, alpha = 69.31 / tMean, 
                                        sigmae = 2 * zSD, t = tMean),
                 g0 = zMax + 2 * zSD),
    
    parPriorMCMC = function(par) {
      dexp(par[1], rate = tMean / 6.931, log = TRUE) +
        dnorm(par[2], zMean, 2 * zSD, log = TRUE) +
        dexp(par[3],  rate = sqrt(tMean / (zVar * 0.6931)), log = TRUE) +
        dnorm(par[4], zMean, 2 * zSD, log = TRUE)
    },
    
    parInitMCMC = function(chainNo, fitML = NULL) {
      if(!is.null(fitML)) {
        parML <- fitML$par
      } else {
        parML <- NULL
      }
      
      init <- rbind(
        c(alpha = 0, theta = zMean, sigma = 1, g0 = zMean),
        parML,
        c(alpha = 0, theta = zMean, sigma = 1, g = zMean)
      )
      
      init[(chainNo - 1) %% nrow(init) + 1, ]
    },
    
    parScaleMCMC = diag(4),
    
    validateSpec = FALSE
  )
  
  if(is.list(spec)) {
    for(name in names(specDefault)) {
      if(is.null(spec[[name]]) & !is.null(specDefault[[name]])) {
        spec[[name]] <- specDefault[[name]]
      } 
    } 
  } else {
    spec <- specDefault
  }
  
  validateSpecPOUMM(spec)
  spec
}

#' @describeIn specPOUMM Fitting a POUMM model with sampling of g0.
#'  Parameter vector is c(alpha, theta, sigma, sigmae, g0).
#' @export
specifyPOUMM_ATSSeG0 <- function(
  z = NULL, tree = NULL,
  zMin = -10, zMean = 0, zMax = 10, zVar = 4, zSD = sqrt(zVar), 
  tMin = 0.1, tMean = 2, tMax = 10, 
  parMapping = NULL, 
  parLower = NULL, parUpper = NULL, 
  g0Prior = NULL,
  parInitML = NULL,
  control = NULL,
  parPriorMCMC = NULL, 
  parInitMCMC = NULL, 
  parScaleMCMC = NULL,
  nSamplesMCMC = 1e5, nAdaptMCMC = nSamplesMCMC, 
  thinMCMC = 100, 
  accRateMCMC = .01, gammaMCMC = 0.50001, nChainsMCMC = 3, 
  samplePriorMCMC = TRUE,
  parallelMCMC = FALSE) {
  
  spec <- list(parMapping = parMapping, 
               parLower = parLower, parUpper = parUpper, 
               g0Prior = g0Prior,
               parInitML = parInitML,
               control = control,
               parPriorMCMC = parPriorMCMC, 
               parInitMCMC = parInitMCMC, 
               parScaleMCMC = parScaleMCMC,
               nSamplesMCMC = nSamplesMCMC, nAdaptMCMC = nAdaptMCMC, 
               thinMCMC = thinMCMC, accRateMCMC = accRateMCMC, 
               gammaMCMC = gammaMCMC, nChainsMCMC = nChainsMCMC, 
               samplePriorMCMC = samplePriorMCMC,
               parallelMCMC = parallelMCMC)
  
  if(validateZTree(z, tree)) {
    zMin <- min(z); zMean <- mean(z); zMax <- max(z); 
    zVar <- var(z); zSD <- sd(z);
    tipTimes <- nodeTimes(tree, tipsOnly = TRUE)
    tMin <- min(tipTimes); tMean <- mean(tipTimes); tMax <- max(tipTimes); 
  }
  
  specDefault <- specifyPOUMM(
    z, tree, zMin, zMean, zMax, zVar, zSD, tMin, tMean, tMax, 
    
    parMapping = function(par) {
      if(is.matrix(par)) {
        colnames(par) <- c("alpha", "theta", "sigma", "sigmae", "g0")
      } else {
        names(par) <- c("alpha", "theta", "sigma", "sigmae", "g0")
      }
      par
    },
    
    parLower = c(alpha = 0, theta = zMin - 2 * (zMax - zMin), 
                 sigma = 0, sigmae = 0, g0 = zMin - 2 * zSD), 
    
    parUpper = c(alpha = 69.31 / tMean, theta = zMax + 2 * (zMax - zMin), 
                 sigma = sigmaOU(H2 = .99, alpha = 69.31 / tMean, 
                                        sigmae = 2 * zSD, t = tMean), 
                 sigmae = 2 * zSD, g0 = zMax + 2 * zSD),
    
    parPriorMCMC = function(par) {
      dexp(par[1], rate = tMean / 6.931, log = TRUE) +
        dnorm(par[2], zMean, 2 * zSD, TRUE) +
        dexp(par[3],  rate = sqrt(tMean / (zVar * 0.6931)), log = TRUE) +
        dexp(par[4], rate = 2 / zSD, log = TRUE) + 
        dnorm(par[5], zMean, 2 * zSD, log = TRUE)
    },
    
    parInitMCMC = function(chainNo, fitML = NULL) {
      if(!is.null(fitML)) {
        parML <- fitML$par
      } else {
        parML <- NULL
      }
      
      init <- rbind(
        c(alpha = 0, theta = zMean, sigma = 1, sigmae = 0, g0 = zMean),
        parML,
        c(alpha = 0, theta = zMean, sigma = 1, sigmae = 1, g0 = zMean)
      )
      
      init[(chainNo - 1) %% nrow(init) + 1, ]
    },
    
    parScaleMCMC = diag(5),
    
    validateSpec = FALSE
  )
  
  if(is.list(spec)) {
    for(name in names(specDefault)) {
      if(is.null(spec[[name]]) & !is.null(specDefault[[name]])) {
        spec[[name]] <- specDefault[[name]]
      } 
    } 
  } else {
    spec <- specDefault
  }
  
  validateSpecPOUMM(spec)
  spec
}


#' @describeIn specPOUMM Specify parameter for fitting a PMM model. 
#'   Parameter vector is c(sigma, sigmae)
#' @export
specifyPMM <- function(
  z = NULL, tree = NULL,
  zMin = -10, zMean = 0, zMax = 10, zVar = 4, zSD = sqrt(zVar), 
  tMin = 0.1, tMean = 2, tMax = 10, 
  parMapping = NULL, 
  parLower = NULL, parUpper = NULL, 
  g0Prior = NULL,
  parInitML = NULL,
  control = NULL,
  parPriorMCMC = NULL, 
  parInitMCMC = NULL, 
  parScaleMCMC = NULL,
  nSamplesMCMC = 1e5, nAdaptMCMC = nSamplesMCMC, 
  thinMCMC = 100, 
  accRateMCMC = .01, gammaMCMC = 0.50001, nChainsMCMC = 3, 
  samplePriorMCMC = TRUE,
  parallelMCMC = FALSE) {

  
  spec <- list(parMapping = parMapping, 
               parLower = parLower, parUpper = parUpper, 
               g0Prior = g0Prior,
               parInitML = parInitML,
               control = control,
               parPriorMCMC = parPriorMCMC, 
               parInitMCMC = parInitMCMC, 
               parScaleMCMC = parScaleMCMC,
               nSamplesMCMC = nSamplesMCMC, nAdaptMCMC = nAdaptMCMC, 
               thinMCMC = thinMCMC, accRateMCMC = accRateMCMC, 
               gammaMCMC = gammaMCMC, nChainsMCMC = nChainsMCMC, 
               samplePriorMCMC = samplePriorMCMC,
               parallelMCMC = parallelMCMC)
  
  if(validateZTree(z, tree)) {
    zMin <- min(z); zMean <- mean(z); zMax <- max(z); 
    zVar <- var(z); zSD <- sd(z);
    tipTimes <- nodeTimes(tree, tipsOnly = TRUE)
    tMin <- min(tipTimes); tMean <- mean(tipTimes); tMax <- max(tipTimes); 
  }
  
  specDefault <- specifyPOUMM(
    z, tree, zMin, zMean, zMax, zVar, zSD, tMin, tMean, tMax, 
    
    parMapping = function(par) {
      if(is.matrix(par)) {
        atsseg0 <- cbind(0, 0, par[, 1:2, drop = FALSE], NA) 
        colnames(atsseg0) <- c("alpha", "theta", "sigma", "sigmae", "g0")
      } else {
        atsseg0 <- c(0, 0, par[1:2], NA) 
        names(atsseg0) <- c("alpha", "theta", "sigma", "sigmae", "g0")
      }
      atsseg0
    },
    
    parLower = c(sigma = 0, sigmae = 0),
    parUpper = c(
      sigma = sigmaOU(H2 = .99, alpha = 0, sigmae = 2 * zSD, t = tMean), 
      sigmae = 2 * zSD),
    
    parPriorMCMC = function(par) {
      dexp(par[1],  rate = sqrt(tMean / (zVar * 0.6931)), log = TRUE) +
        dexp(par[2], rate = 2 / zSD, log = TRUE)
    },

    parInitMCMC = function(chainNo, fitML = NULL) {
      if(!is.null(fitML)) {
        parML <- fitML$par
      } else {
        parML <- NULL
      }
      
      init <- rbind(
        c(sigma = 1, sigmae = 0),
        parML,
        c(sigma = 1, sigmae = 1)
      )
      
      init[(chainNo - 1) %% nrow(init) + 1, ]
    },
    parScaleMCMC = diag(2), 
    
    validateSpec = FALSE
  )
  
  if(is.list(spec)) {
    for(name in names(specDefault)) {
      if(is.null(spec[[name]]) & !is.null(specDefault[[name]])) {
        spec[[name]] <- specDefault[[name]]
      } 
    } 
  } else {
    spec <- specDefault
  }
  
  validateSpecPOUMM(spec)
  spec
}

#' @describeIn specPOUMM Specify parameter for fitting a PMM model with
#'  sampling of g0. Parameter vector is c(sigma, sigmae, g0).
#' @export
specifyPMM_SSeG0 <- function(
  z = NULL, tree = NULL,
  zMin = -10, zMean = 0, zMax = 10, zVar = 4, zSD = sqrt(zVar), 
  tMin = 0.1, tMean = 2, tMax = 10, 
  parMapping = NULL, 
  parLower = NULL, parUpper = NULL, 
  g0Prior = NULL,
  parInitML = NULL,
  control = NULL,
  parPriorMCMC = NULL, 
  parInitMCMC = NULL, 
  parScaleMCMC = NULL,
  nSamplesMCMC = 1e5, nAdaptMCMC = nSamplesMCMC, 
  thinMCMC = 100, 
  accRateMCMC = .01, gammaMCMC = 0.50001, nChainsMCMC = 3, 
  samplePriorMCMC = TRUE,
  parallelMCMC = FALSE) {
  
  spec <- list(parMapping = parMapping, 
               parLower = parLower, parUpper = parUpper, 
               g0Prior = g0Prior,
               parInitML = parInitML,
               control = control,
               parPriorMCMC = parPriorMCMC, 
               parInitMCMC = parInitMCMC, 
               parScaleMCMC = parScaleMCMC,
               nSamplesMCMC = nSamplesMCMC, nAdaptMCMC = nAdaptMCMC, 
               thinMCMC = thinMCMC, accRateMCMC = accRateMCMC, 
               gammaMCMC = gammaMCMC, nChainsMCMC = nChainsMCMC, 
               samplePriorMCMC = samplePriorMCMC,
               parallelMCMC = parallelMCMC)
  
  if(validateZTree(z, tree)) {
    zMin <- min(z); zMean <- mean(z); zMax <- max(z); 
    zVar <- var(z); zSD <- sd(z);
    tipTimes <- nodeTimes(tree, tipsOnly = TRUE)
    tMin <- min(tipTimes); tMean <- mean(tipTimes); tMax <- max(tipTimes); 
  }
  
  specDefault <- specifyPOUMM(
    z, tree, zMin, zMean, zMax, zVar, zSD, tMin, tMean, tMax, 
    
    parMapping = function(par) {
      if(is.matrix(par)) {
        atsseg0 <- cbind(0, 0, par[, 1:3, drop = FALSE]) 
        colnames(atsseg0) <- c("alpha", "theta", "sigma", "sigmae", "g0")
      } else {
        atsseg0 <- c(0, 0, par[1:3]) 
        names(atsseg0) <- c("alpha", "theta", "sigma", "sigmae", "g0")
      }
      atsseg0
    },
    
    parLower = c(sigma = 0, sigmae = 0, g0 = zMin - 2 * zSD),
    parUpper = c(
      sigma = sigmaOU(H2 = .99, alpha = 0, sigmae = 2 * zSD, t = tMean), 
      sigmae = 2 * zSD, g0 = zMax + 2 * zSD),
    
    parPriorMCMC = function(par) {
      dexp(par[1],  rate = sqrt(tMean / (zVar * 0.6931)), log = TRUE) +
        dexp(par[2], rate = 2 / zSD, log = TRUE) + 
        dnorm(par[3], zMean, 2 * zSD, log = TRUE)
    },
    
    parInitMCMC = function(chainNo, fitML = NULL) {
      if(!is.null(fitML)) {
        parML <- fitML$par
      } else {
        parML <- NULL
      }
      
      init <- rbind(
        c(sigma = 1, sigmae = 0, g0 = zMean),
        parML,
        c(sigma = 1, sigmae = 1, g0 = zMean)
      )
      
      init[(chainNo - 1) %% nrow(init) + 1, ]
    },
    parScaleMCMC = diag(3), 
    
    validateSpec = FALSE
  )
  
  if(is.list(spec)) {
    for(name in names(specDefault)) {
      if(is.null(spec[[name]]) & !is.null(specDefault[[name]])) {
        spec[[name]] <- specDefault[[name]]
      } 
    } 
  } else {
    spec <- specDefault
  }
  
  validateSpecPOUMM(spec)
  spec
}

#' @describeIn specPOUMM Fitting a POUMM model with a uniform prior for
#'  the phylogenetic heritability at mean root-tip distance. Parameter vector is
#'  c(alpha, theta, H2tMean, sigmae).
#' @export
specifyPOUMM_ATH2tMeanSe <- function(
  z = NULL, tree = NULL,
  zMin = -10, zMean = 0, zMax = 10, zVar = 4, zSD = sqrt(zVar), 
  tMin = 0.1, tMean = 2, tMax = 10, 
  parMapping = NULL, 
  parLower = NULL, parUpper = NULL, 
  g0Prior = NULL,
  parInitML = NULL,
  control = NULL,
  parPriorMCMC = NULL, 
  parInitMCMC = NULL, 
  parScaleMCMC = NULL,
  nSamplesMCMC = 1e5, nAdaptMCMC = nSamplesMCMC, 
  thinMCMC = 100, 
  accRateMCMC = .01, gammaMCMC = 0.50001, nChainsMCMC = 3, 
  samplePriorMCMC = TRUE,
  parallelMCMC = FALSE) {
  
  spec <- list(parMapping = parMapping, 
               parLower = parLower, parUpper = parUpper, 
               g0Prior = g0Prior,
               parInitML = parInitML,
               control = control,
               parPriorMCMC = parPriorMCMC, 
               parInitMCMC = parInitMCMC, 
               parScaleMCMC = parScaleMCMC,
               nSamplesMCMC = nSamplesMCMC, nAdaptMCMC = nAdaptMCMC, 
               thinMCMC = thinMCMC, accRateMCMC = accRateMCMC, 
               gammaMCMC = gammaMCMC, nChainsMCMC = nChainsMCMC, 
               samplePriorMCMC = samplePriorMCMC,
               parallelMCMC = parallelMCMC)
  
  if(validateZTree(z, tree)) {
    zMin <- min(z); zMean <- mean(z); zMax <- max(z); 
    zVar <- var(z); zSD <- sd(z);
    tipTimes <- nodeTimes(tree, tipsOnly = TRUE)
    tMin <- min(tipTimes); tMean <- mean(tipTimes); tMax <- max(tipTimes); 
  }
  
  specDefault <- specifyPOUMM(
    z, tree, zMin, zMean, zMax, zVar, zSD, tMin, tMean, tMax, 
    
    parMapping = function(par) {
      if(is.matrix(par)) {
        par[, 3] <- sigmaOU(par[, 3], par[, 1], par[, 4], tMean)
        par <- cbind(par, NA)
        
        colnames(par) <- c("alpha", "theta", "sigma", "sigmae", "g0")
      } else {
        par[3] <- sigmaOU(par[3], par[1], par[4], tMean)
        par <- c(par, NA)
        
        names(par) <- c("alpha", "theta", "sigma", "sigmae", "g0")
      }
      par
    },
    
    parLower = c(alpha = 0, theta = zMin - 2 * (zMax - zMin), 
                 H2tMean = 0, sigmae = 0), 
    
    parUpper = c(alpha = 69.31 / tMean, theta = zMax + 2 * (zMax - zMin), 
                 H2tMean = .99, sigmae = 2 * zSD),
  
    parPriorMCMC = function(par) {
      dexp(par[1], rate = tMean / 6.931, log = TRUE) +
        dnorm(par[2], zMean, 2 * zSD, log = TRUE) +
        dunif(par[3], min = 0, max = 1, log = TRUE) +
        dexp(par[4], rate = 2 / zSD, log = TRUE)
    },
    
    parInitMCMC = function(chainNo, fitML = NULL) {
      if(!is.null(fitML)) {
        parML <- fitML$par
      } else {
        parML <- NULL
      }
      
      init <- rbind(
        c(alpha = 0, theta = zMean, H2tMean = .9, sigmae = 0),
        parML,
        c(alpha = 0, theta = zMean, H2tMean = .1, sigmae = 1)
      )
      
      init[(chainNo - 1) %% nrow(init) + 1, ]
    },
    
    parScaleMCMC = diag(4), 
    
    validateSpec = FALSE
  )
  
  if(is.list(spec)) {
    for(name in names(specDefault)) {
      if(is.null(spec[[name]]) & !is.null(specDefault[[name]])) {
        spec[[name]] <- specDefault[[name]]
      } 
    } 
  } else {
    spec <- specDefault
  }
  
  validateSpecPOUMM(spec)
  spec
}



#' @describeIn specPOUMM Fitting a POUMM model with a uniform prior for
#'  the phylogenetic heritability at mean root-tip with sampling of g0.
#'  Parameter vector is c(alpha, theta, H2tMean, sigmae, g0).
#' @export
specifyPOUMM_ATH2tMeanSeG0 <- function(
  z = NULL, tree = NULL,
  zMin = -10, zMean = 0, zMax = 10, zVar = 4, zSD = sqrt(zVar), 
  tMin = 0.1, tMean = 2, tMax = 10, 
  parMapping = NULL, 
  parLower = NULL, parUpper = NULL, 
  g0Prior = NULL,
  parInitML = NULL,
  control = NULL,
  parPriorMCMC = NULL, 
  parInitMCMC = NULL, 
  parScaleMCMC = NULL,
  nSamplesMCMC = 1e5, nAdaptMCMC = nSamplesMCMC, 
  thinMCMC = 100, 
  accRateMCMC = .01, gammaMCMC = 0.50001, nChainsMCMC = 3, 
  samplePriorMCMC = TRUE,
  parallelMCMC = FALSE) {
  
  spec <- list(parMapping = parMapping, 
               parLower = parLower, parUpper = parUpper, 
               g0Prior = g0Prior,
               parInitML = parInitML,
               control = control,
               parPriorMCMC = parPriorMCMC, 
               parInitMCMC = parInitMCMC, 
               parScaleMCMC = parScaleMCMC,
               nSamplesMCMC = nSamplesMCMC, nAdaptMCMC = nAdaptMCMC, 
               thinMCMC = thinMCMC, accRateMCMC = accRateMCMC, 
               gammaMCMC = gammaMCMC, nChainsMCMC = nChainsMCMC, 
               samplePriorMCMC = samplePriorMCMC,
               parallelMCMC = parallelMCMC)
  
  if(validateZTree(z, tree)) {
    zMin <- min(z); zMean <- mean(z); zMax <- max(z); 
    zVar <- var(z); zSD <- sd(z);
    tipTimes <- nodeTimes(tree, tipsOnly = TRUE)
    tMin <- min(tipTimes); tMean <- mean(tipTimes); tMax <- max(tipTimes); 
  }
  
  specDefault <- specifyPOUMM(
    z, tree, zMin, zMean, zMax, zVar, zSD, tMin, tMean, tMax, 
    
    parMapping = function(par) {
      if(is.matrix(par)) {
        par[, 3] <- sigmaOU(par[, 3], par[, 1], par[, 4], tMean)
        
        colnames(par) <- c("alpha", "theta", "sigma", "sigmae", "g0")
      } else {
        par[3] <- sigmaOU(par[3], par[1], par[4], tMean)
        
        names(par) <- c("alpha", "theta", "sigma", "sigmae", "g0")
      }
      par
    },
    
    parLower = c(alpha = 0, theta = zMin - 2 * (zMax - zMin), 
                 H2tMean = 0, sigmae = 0, g0 = zMin - 2 * zSD), 
    
    parUpper = c(alpha = 69.31 / tMean, theta = zMax + 2 * (zMax - zMin), 
                 H2tMean = .99, sigmae = 2 * zSD, g0 = zMax + 2 * zSD),
    
    parPriorMCMC = function(par) {
      dexp(par[1], rate = tMean / 6.931, log = TRUE) +
        dnorm(par[2], zMean, 2 * zSD, TRUE) +
        dunif(par[3], min = 0, max = 1, log = TRUE) +
        dexp(par[4], rate = 2 / zSD, log = TRUE) + 
        dnorm(par[5], zMean, 2 * zSD, log = TRUE)
    },
    
    parInitMCMC = function(chainNo, fitML = NULL) {
      if(!is.null(fitML)) {
        parML <- fitML$par
      } else {
        parML <- NULL
      }
      
      init <- rbind(
        c(alpha = 0, theta = zMean, H2tMean = .9, sigmae = 0, g0 = zMean),
        parML,
        c(alpha = 0, theta = zMean, H2tMean = .1, sigmae = 1, g0 = zMean)
      )
      
      init[(chainNo - 1) %% nrow(init) + 1, ]
    },
    
    parScaleMCMC = diag(5), 
    
    validateSpec = FALSE
  )
  
  if(is.list(spec)) {
    for(name in names(specDefault)) {
      if(is.null(spec[[name]]) & !is.null(specDefault[[name]])) {
        spec[[name]] <- specDefault[[name]]
      } 
    } 
  } else {
    spec <- specDefault
  }
  
  validateSpecPOUMM(spec)
  spec
}



#' @describeIn specPOUMM Fitting a PMM model with a uniform prior for
#'  the phylogenetic heritability at mean root-tip distance. Parameter vector is
#'  c(H2tMean, sigmae).
#' @export
specifyPMM_H2tMeanSe <- function(
  z = NULL, tree = NULL,
  zMin = -10, zMean = 0, zMax = 10, zVar = 4, zSD = sqrt(zVar), 
  tMin = 0.1, tMean = 2, tMax = 10, 
  parMapping = NULL, 
  parLower = NULL, parUpper = NULL, 
  g0Prior = NULL,
  parInitML = NULL,
  control = NULL,
  parPriorMCMC = NULL, 
  parInitMCMC = NULL, 
  parScaleMCMC = NULL,
  nSamplesMCMC = 1e5, nAdaptMCMC = nSamplesMCMC, 
  thinMCMC = 100, 
  accRateMCMC = .01, gammaMCMC = 0.50001, nChainsMCMC = 3, 
  samplePriorMCMC = TRUE,
  parallelMCMC = FALSE) {
  
  spec <- list(parMapping = parMapping, 
               parLower = parLower, parUpper = parUpper, 
               g0Prior = g0Prior,
               parInitML = parInitML,
               control = control,
               parPriorMCMC = parPriorMCMC, 
               parInitMCMC = parInitMCMC, 
               parScaleMCMC = parScaleMCMC,
               nSamplesMCMC = nSamplesMCMC, nAdaptMCMC = nAdaptMCMC, 
               thinMCMC = thinMCMC, accRateMCMC = accRateMCMC, 
               gammaMCMC = gammaMCMC, nChainsMCMC = nChainsMCMC, 
               samplePriorMCMC = samplePriorMCMC,
               parallelMCMC = parallelMCMC)
  
  if(validateZTree(z, tree)) {
    zMin <- min(z); zMean <- mean(z); zMax <- max(z); 
    zVar <- var(z); zSD <- sd(z);
    tipTimes <- nodeTimes(tree, tipsOnly = TRUE)
    tMin <- min(tipTimes); tMean <- mean(tipTimes); tMax <- max(tipTimes); 
  }
  
  specDefault <- specifyPOUMM(
    z, tree, zMin, zMean, zMax, zVar, zSD, tMin, tMean, tMax, 
    
    parMapping = function(par) {
      if(is.matrix(par)) {
        par <- cbind(0, 0, par[, 1:2, drop = FALSE], NA)
        par[, 3] <- sigmaOU(par[, 3], 0, par[, 4], tMean)
        
        colnames(par) <- c("alpha", "theta", "sigma", "sigmae", "g0")
      } else {
        par <- c(0, 0, par, NA)
        par[3] <- sigmaOU(par[3], 0, par[4], tMean)
        
        names(par) <- c("alpha", "theta", "sigma", "sigmae", "g0")
      }
      par
    },
    
    parLower = c(H2tMean = 0, sigmae = 0), 
    
    parUpper = c(H2tMean = .99, sigmae = 2 * zSD),
    
    parPriorMCMC = function(par) {
      dunif(par[1], min = 0, max = 1, log = TRUE) +
        dexp(par[2], rate = 2 / zSD, log = TRUE)
    },
    
    parInitMCMC = function(chainNo, fitML = NULL) {
      if(!is.null(fitML)) {
        parML <- fitML$par
      } else {
        parML <- NULL
      }
      
      init <- rbind(
        c(H2tMean = .9, sigmae = 0),
        parML,
        c(H2tMean = .1, sigmae = 1)
      )
      
      init[(chainNo - 1) %% nrow(init) + 1, ]
    },
    
    parScaleMCMC = diag(2), 
    
    validateSpec = FALSE
  )
  
  if(is.list(spec)) {
    for(name in names(specDefault)) {
      if(is.null(spec[[name]]) & !is.null(specDefault[[name]])) {
        spec[[name]] <- specDefault[[name]]
      } 
    } 
  } else {
    spec <- specDefault
  }
  
  validateSpecPOUMM(spec)
  spec
}


#' @describeIn specPOUMM Fitting a PMM model with a uniform prior for
#'  the phylogenetic heritability at mean root-tip distance with sampling of G0.
#'  Parameter vector is c(H2tMean, sigmae, g0).
#' @export
specifyPMM_H2tMeanSeG0 <- function(
  z = NULL, tree = NULL,
  zMin = -10, zMean = 0, zMax = 10, zVar = 4, zSD = sqrt(zVar), 
  tMin = 0.1, tMean = 2, tMax = 10, 
  parMapping = NULL, 
  parLower = NULL, parUpper = NULL, 
  g0Prior = NULL,
  parInitML = NULL,
  control = NULL,
  parPriorMCMC = NULL, 
  parInitMCMC = NULL, 
  parScaleMCMC = NULL,
  nSamplesMCMC = 1e5, nAdaptMCMC = nSamplesMCMC, 
  thinMCMC = 100, 
  accRateMCMC = .01, gammaMCMC = 0.50001, nChainsMCMC = 3, 
  samplePriorMCMC = TRUE,
  parallelMCMC = FALSE) {
  
  spec <- list(parMapping = parMapping, 
               parLower = parLower, parUpper = parUpper, 
               g0Prior = g0Prior,
               parInitML = parInitML,
               control = control,
               parPriorMCMC = parPriorMCMC, 
               parInitMCMC = parInitMCMC, 
               parScaleMCMC = parScaleMCMC,
               nSamplesMCMC = nSamplesMCMC, nAdaptMCMC = nAdaptMCMC, 
               thinMCMC = thinMCMC, accRateMCMC = accRateMCMC, 
               gammaMCMC = gammaMCMC, nChainsMCMC = nChainsMCMC, 
               samplePriorMCMC = samplePriorMCMC,
               parallelMCMC = parallelMCMC)
  
  if(validateZTree(z, tree)) {
    zMin <- min(z); zMean <- mean(z); zMax <- max(z); 
    zVar <- var(z); zSD <- sd(z);
    tipTimes <- nodeTimes(tree, tipsOnly = TRUE)
    tMin <- min(tipTimes); tMean <- mean(tipTimes); tMax <- max(tipTimes); 
  }
  
  specDefault <- specifyPOUMM(
    z, tree, zMin, zMean, zMax, zVar, zSD, tMin, tMean, tMax, 
    
    parMapping = function(par) {
      if(is.matrix(par)) {
        par <- cbind(0, 0, par)
        par[, 3] <- sigmaOU(par[, 3], 0, par[, 4], tMean)
        
        colnames(par) <- c("alpha", "theta", "sigma", "sigmae", "g0")
      } else {
        par <- c(0, 0, par)
        par[3] <- sigmaOU(par[3], par[1], par[4], tMean)
        
        names(par) <- c("alpha", "theta", "sigma", "sigmae", "g0")
      }
      par
    },
    
    parLower = c(H2tMean = 0, sigmae = 0, g0 = zMin - 2 * zSD), 
    
    parUpper = c(H2tMean = .99, sigmae = 2 * zSD, g0 = zMax + 2 * zSD),
    
    parPriorMCMC = function(par) {
      dunif(par[1], min = 0, max = 1, log = TRUE) +
        dexp(par[2], rate = 2 / zSD, log = TRUE) + 
        dnorm(par[3], zMean, 2 * zSD, log = TRUE)
      
    },
    
    parInitMCMC = function(chainNo, fitML = NULL) {
      if(!is.null(fitML)) {
        parML <- fitML$par
      } else {
        parML <- NULL
      }
      
      init <- rbind(
        c(H2tMean = .9, sigmae = 0, g0 = zMean),
        parML,
        c(H2tMean = .1, sigmae = 1, g0 = zMean)
      )
      
      init[(chainNo - 1) %% nrow(init) + 1,]
    },
    
    parScaleMCMC = diag(3), 
    
    validateSpec = FALSE
  )
  
  if(is.list(spec)) {
    for(name in names(specDefault)) {
      if(is.null(spec[[name]]) & !is.null(specDefault[[name]])) {
        spec[[name]] <- specDefault[[name]]
      } 
    } 
  } else {
    spec <- specDefault
  }
  
  validateSpecPOUMM(spec)
  spec
}


######### Validate specification ###########
#' Validate phenotypic values and phylogenetic tree
#' @param z trait (phenotypic) values at the tips of the tree
#' @param tree A phylo object with the same number of tips as the length of z.
#' @return The function either returns TRUE or exits with an error message if it
#'  finds a problem with the specificaiton.
validateZTree <- function(z, tree) {
  if(!is.null(z) & !is.null(tree)) {
    if(!is.vector(z, mode="double")) {
      stop("The trait-vector z should be a vector of mode 'double'!")
    }
    if(length(z) == 0) {
      stop("The trait-vector z is empty, but should contain at least 1 element!")
    }
    if(any(is.na(z)) | any(is.infinite(z))) {
      stop("The trait vector z contains infinite or NA values. All trait-values should be finite!")
    }
    if(!("phylo" %in% class(tree))) {
      stop("tree must be a phylo object!")
    }
    if(length(tree$tip.label) != length(z)) {
      stop("z should have the same length as the number of tips in tree!")
    } 
    if(any(tree$edge.length <= 0)) {
      stop('All edge lengths in tree should be positive!')
    }
    if(!is.null(names(z))) {
      if(!all(names(z) == tree$tip.label)) {
        stop("Some of the names in the trait-vector do not correspond to tip-labels. names(z) should be either NULL or it should be identical (in the same order) as tree$tip.label.")  
      }
    }
    TRUE  
  } else {
    FALSE
  }
}
#' Validate a POUMM specification
#' @param spec A list object returned by one of the specifyPOUMM or specifyPMM
#'  functions with possibly modified entries afterwards.
#' 
#' @return The function either returns TRUE or exits with an error message if it
#'  finds a problem with the specificaiton.
#'  
#' @export
validateSpecPOUMM <- function(spec) {
  with(spec, {
    if(nChainsMCMC <= 0) {
      stop(paste0("nChainsMCMC should be greater or equal to 1 but was ", 
                  nChainsMCMC, 
                  ". To disable MCMC-fit, set doMCMC = FALSE in the POUMM call."))
    }
    if(is.numeric(nAdaptMCMC) & nAdaptMCMC < 0) {
      stop("nAdaptMCMC should be logical or a non-negative integer.")
    }
    if(nSamplesMCMC < 0) {
      stop("nSamplesMCMC should be a non-negative integer (0 disables MCMC run).")
    }
    if(thinMCMC <= 0) {
      stop("thinMCMC should be a positive integer (1 means no skipping).")
    }
    if(accRateMCMC <= 0 | accRateMCMC > 1) {
      stop("accRateMCMC should be a positive number between 0 and 1.")
    } 
    if(gammaMCMC <= 0.5 | gammaMCMC > 1) {
      stop("gammaMCMC should be in (0.5,1].")
    }
    if(nChainsMCMC == 1 & samplePriorMCMC) {
      warning("Only one MCMC to be run which samples from the prior.")
    }
    
    
    listParInitMCMC <- lapply(1:spec$nChainsMCMC, function(i) parInitMCMC(i))
    parNames <- names(listParInitMCMC[[1]])
    parLen <- length(listParInitMCMC[[1]])
    
    if(!is.null(parInitML)) {
      if(is.list(parInitML)) {
        sapply(parInitML, function(v) {
          if(!is.vector(v) | length(v) != parLen | 
             !identical(names(v), parNames)) {
            stop("If passing a list for parInitML, this should contain vectors of the same length and with the same names as each vector returned by parInitMCMC.")
          }
        } )
      } else if(!is.vector(parInitML) | length(parInitML) != parLen | 
         !identical(names(parInitML), parNames)) {
        stop("parInitML should be a vector of the same length and with the same names as each vector returned by parInitMCMC.")
      }
    }
    
    if(is.null(parLower) | is.null(parUpper) | 
       any(is.na(parLower)) | any(is.na(parUpper))) {
      cat("parLower:")
      print(parLower)
      cat("parUpper:")
      print(parUpper)
      
      
      stop("NULL parLower and/or parUpper or some NA or NaNs found in parLower
           and/or parUpper. Did you forget
           specifying some of the parameters zMin, zMean, zMax, zVar, tMin,
           tMean, tMax? If yes, you can fix them or specify parLower and
           parUpper explicitly.")
    }
    if(length(parLower) != parLen | length(parUpper) != parLen | 
       !identical(parNames, names(parLower)) | !identical(parNames, names(parUpper))) {
      warning("parInitMCMC returns vectors of different length or with different names from parLower and parUpper.")
    } 
    
    for(i in 1:spec$nChainsMCMC) {
      parInitMCMC <- listParInitMCMC[[i]]
      if(is.null(parInitMCMC) | !is.vector(parInitMCMC) | length(parInitMCMC) != parLen |
         is.null(names(parInitMCMC)) | !identical(names(parInitMCMC), parNames)) {
        print(parInitMCMC)
        stop(paste0("parInitMCMC should return a named vector with the same length and names for each chain index.", " Check parInitMCMC(",i,")."))
      }
      
      atsseg0 <- parMapping(parInitMCMC)
      if(!is.vector(atsseg0) | is.null(names(atsseg0)) | length(atsseg0) != 5 |
         !identical(names(atsseg0), c("alpha", "theta", "sigma", "sigmae", "g0")) |
         atsseg0[1] < 0 | atsseg0[3] < 0 | atsseg0[4] < 0 | any(is.na(atsseg0[1:4]))) {
        print("parInitMCMC : ")
        print(parInitMCMC)
        print(" mapped to : ")
        print(atsseg0)
        stop(paste0("When given a vector, parMapping should return a vector ",
                    "with named elements alpha>=0, theta, sigma>=0, sigmae>=0 and g0.", 
                    "Only g0 can be NA or NaN which has a specical meaning (described in doc for parMapping). Check parInitMCMC(",i,")."))
      }
      
      prior <- try(parPriorMCMC(parInitMCMC), silent = TRUE)
      if( (!"numeric" %in% class(prior)) | is.na(prior) ) {
        print("parInitMCMC : ")
        print(parInitMCMC)
        stop(paste0("parPriorMCMC should return a finite number but returned:", 
                    toString(prior)))
      }
      
      parInitMat <- matrix(parInitMCMC, nrow = 1)
      colnames(parInitMat) <- names(parInitMCMC)
      atsseg0Mat <- parMapping(parInitMat)
      
      
      if(!is.matrix(atsseg0Mat) | is.null(colnames(atsseg0Mat)) | 
         length(atsseg0Mat) != 5 |
         !identical(colnames(atsseg0Mat), c("alpha", "theta", "sigma", "sigmae", "g0")) |
         any(atsseg0Mat[, 1] < 0) | any(atsseg0Mat[, 3] < 0) | 
         any(atsseg0Mat[, 4] < 0) | any(is.na(atsseg0Mat[, 1:4]))) {
        
        cat("parInitMCMC : \n")
        print(parInitMat)
        cat(" mapped to : \n")
        print(atsseg0Mat)
        stop(paste0("When given a matrix, parMapping should return a matrix ",
                    "with named columns alpha>=0, theta, sigma>=0, sigmae>=0 and g0.", 
                    "Only g0 can be NA or NaN which has a specical meaning ",
                    "(described in doc for parMapping)."))
      }
    }
  })
  TRUE
}
