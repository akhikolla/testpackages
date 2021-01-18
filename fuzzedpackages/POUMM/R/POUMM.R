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


# Implementation of the POUMM likelihood and heritability estimators

#'@title The Phylogenetic (Ornstein-Uhlenbeck) Mixed Model
#'  
#'@description This is the high-level entry point to the POUMM method. The POUMM
#'  function fits the POUMM method to a tree and observed trait-values at its tips
#'   and returns an object of class "POUMM".
#'  
#'@param z Either a numeric vector containing the phenotypic values at the tips 
#'  of tree or a named list containing named elements z - a numeric vector and 
#'  tree - a phylo object (it is possible to specify different element names 
#'  using the arguments zName and treeName).
#'@param tree A phylo object or NULL in case z is a list.
#'@param se A non-negative numerical vector (or single number) indicating known 
#'measurement standard error (defaults to 0). Note the elements of this vector 
#'are assumed to describe the measurement error at individual nodes independent
#'of the environmental contribution (described by the parameter sigmae). The total
#'error standard deviation is thus sqrt(sigmae2+se^2).
#'@param zName,treeName Character strings used when the parameter z is a list; 
#'  indicate the names in the list of the values-vector and the tree. Default: 
#'  'z' and 'tree'.
#'@param parDigits Integer specifying rounding to be done on the parameter 
#'  vector before likelihood calculation. Defaults to 6 decimal digits. This can
#'  be useful during maximum likelihood optimization to prevent likelihood 
#'  calculation on very small but positive values of alpha, but should be used 
#'  with caution since specifying a small number of digits, i.e. 2 or 3 can
#'  result in an infinite loop during optim. Specify a negative number
#'  to disable rounding.
#'  
#'@param usempfr integer indicating if and how mpfr should be used for small 
#'  parameter values (`any(c(alpha, sigma, sigmae) < 0.01)`). Using the mpfr 
#'  package can be forced by specifying an integer greater or equal to 2. 
#'  Setting usempfr=0 (default) causes high precision likelihood 
#'  calculation to be done on each encounter of parameters with at least 1 bigger
#'  log-likelihood value than any of the currently found
#'  maximum log-likelihood or the previously calculated log-likelihood value
#'  Requires the Rmpfr package. Note that using mpfr may increase the time for 
#'  one likelihood calculation more than 100-fold. Set usempfr to -1 or less
#'  to completely disable Rmpfr functionality. 
#'  
#'@param useCpp Logical indicating whether C++ likelihood calculation should be 
#'  used for faster vector operations. Defaults to TRUE. Since the C++ likelihood
#'  implementation does not support mpfr, useCpp gets disabled when usempfr is 
#'  bigger than 0.
#'
#'@param ... additional arguments passed to the `likPOUMMGivenTreeVTips()` function 
#'  (`?dVGivenTreeOU` for details).
#'  
#'@param spec A named list specifying how the ML and MCMC fit should be done. 
#'  See `?specifyPOUMM`.
#'  
#'@param doMCMC Deprecated - replaced by specifying nSamplesMCMC as a member
#'  of spec instead (see `?specifyPOUMM`). 
#'  logical: should a MCMC fit be performed. An MCMC fit provides a 
#'  sample from the posterior distribution of the parameters given a prior 
#'  distribution and the data. Unlike the ML-fit, it allows to estimate 
#'  confidence intervals for the estimated parameters. This argument is TRUE by 
#'  default. The current implementation uses the adaptive 
#'  Metropolis sampler from the package `adaptMCMC` written by Andreas 
#'  Scheidegger. To obtain meaningful estimates MCMC may need to run for several
#'  millions of iterations (parameter nSamplesMCMC set to 1e5 by default). See 
#'  parameters ending at MCMC in `?specifyPOUMM` for details.
#'
#'@param likPOUMM_lowLevelFun the low-level function used for POUMM - likelihood 
#'calculation. Default value is likPOUMMGivenTreeVTipsC.
#'
#'@param verbose,debug Logical flags indicating whether to print informative 
#'  and/or debug information on the standard output (both are set to to FALSE by
#'  default).
#'  
#'@return An object of S3 class 'POUMM'. This object can be analyzed using 
#'  S3 generic functions: \code{\link{summary}}, 
#'  \code{\link{plot}}, \code{\link{AIC}}, \code{\link{BIC}}, \code{\link{coef}},
#'  \code{\link{logLik}}, \code{\link{fitted}}.
#' 
#' @seealso \code{\link{specifyPOUMM}} for parametrizations and custom settings
#'  of the POUMM fit.
#'  
#' @examples 
#' \dontrun{
#' # Please, read the package vignette for more detailed examples.
#' N <- 500
#' tr <- ape::rtree(N)
#' z <- rVNodesGivenTreePOUMM(tr, 0, 2, 3, 1, 1)[1:N]
#' fit <- POUMM(z, tr, spec = specifyPOUMM(nSamplesMCMC = 5e4))
#' plot(fit)
#' summary(fit)
#' AIC(fit)
#' BIC(fit)
#' coef(fit)
#' logLik(fit)
#' fitted(fit)
#' plot(resid(fit))
#' abline(h=0)
#' 
#' # fit PMM to the same data and do a likelihood ratio test
#' fitPMM <- POUMM(z, tr, spec = specifyPMM(nSamplesMCMC = 5e4))
#' lmtest::lrtest(fitPMM, fit)
#' }
#' 
#'@references 
#'  Mitov, V., and Stadler, T. (2017). Fast and Robust Inference of Phylogenetic Ornstein-Uhlenbeck Models Using Parallel Likelihood Calculation. bioRxiv, 115089. 
#'  https://doi.org/10.1101/115089
#'  Mitov, V., & Stadler, T. (2017). Fast Bayesian Inference of Phylogenetic Models Using Parallel Likelihood Calculation and Adaptive Metropolis Sampling. Systematic Biology, 235739. http://doi.org/10.1101/235739  
#'  Vihola, M. (2012). Robust adaptive Metropolis algorithm with coerced 
#'  acceptance rate. Statistics and Computing, 22(5), 997-1008. 
#'  http://doi.org/10.1007/s11222-011-9269-5 
#'  
#'  Scheidegger, A. (2012). adaptMCMC: Implementation of a generic adaptive 
#'  Monte Carlo Markov Chain sampler. 
#'  http://CRAN.R-project.org/package=adaptMCMC
#'
#'@importFrom stats var sd rnorm dnorm dexp rexp dunif runif nobs
#'@useDynLib POUMM
#' 
#' @export
#'
POUMM <- function(
  z, tree, se = 0, zName = 'z', treeName = 'tree', 
  parDigits = 6, usempfr = 0, useCpp = TRUE,
  ..., 
  spec = NULL, doMCMC = TRUE,
  likPOUMM_lowLevelFun = likPOUMMGivenTreeVTipsC,
  verbose = FALSE, debug=FALSE) {
  
  ###### Verify input data ######
  if(is.list(z)) {
    p <- z
    z <- p[[zName]]
    if(is.null(z)) {
      z <- p[['v']]
    }
    tree <- p[[treeName]]
    
    if(is.null(z) | is.null(tree)) {
      stop('If a list is supplied as argument z, this list should contain a ',
           'vector of trait values named "z" or zName and a phylo-object named ',
           '"tree" or treeName')
    }
    if(any(is.na(z)) | any(is.infinite(z))) {
      stop('Check z for infinite or NA values!')
    }
    if(any(tree$edge.length <= 0) | any(is.infinite(tree$edge.length)) | 
       any(is.na(tree$edge.length))) {
      stop('Check the tree for non-finite or non-positive edge-lengths!')
    }
  }
  
  ######## Caching pruneInfo for faster likelihood calculations
  if(!validateZTree(z, tree)) {
    stop("Invalid z and/or tree.")
  }
  
  pruneInfo <- pruneTree(tree, z, se)
  
  tTips <- nodeTimes(tree, tipsOnly = TRUE)
  
  
  ######## Default POUMM spec ###########
  if(is.function(spec)) {
    spec <- do.call(spec, 
                    list(z = z, tree = tree, 
                         zMin = min(z), zMean = mean(z), zMax = max(z), 
                         zVar = var(z), zSD = sd(z), 
                         tMin = min(tTips), tMean = mean(tTips), tMax = max(tTips)))
  } else {
    spec <- 
      do.call(specifyPOUMM, 
              c(list(z = z, tree = tree, 
                     zMin = min(z), zMean = mean(z), zMax = max(z), 
                     zVar = var(z), tMin = min(tTips), tMean = mean(tTips),
                     tMax = max(tTips)), spec))
  }
  
  g0Lower <- spec$parMapping(spec$parLower)['g0']
  
  dof = if(is.na(g0Lower) & !is.nan(g0Lower)) {
    length(spec$parLower) + 1 
  } else {
    length(spec$parLower) 
  }
  
  result <- list(pruneInfo = pruneInfo, 
                 N = length(tree$tip.label), dof = dof, 
                 tMax = max(nodeTimes(tree, tipsOnly = TRUE)),
                 tMean = mean(nodeTimes(tree, tipsOnly = TRUE)), 
                 spec = spec, 
                 ...)
  
  # define a loglik function to be called during likelihood maximization as well
  # as MCMC-sampling. 
  # The argument par is a named numeric vector. This vector is mapped to the 
  # POUMM parameters alpha, theta, sigma, sigmae and g0 using the function 
  # spec$parMapping. If g0 = NA the loglik
  # funciton finds the value of g0 that maximizes the likelihood, i.e. 
  # p(z | tree, alpha, theta, sigma, sigmae, g0), or if g0Prior specifiies a 
  # normal distribution with mean g0Prior$mean and variance g0Prior$var, 
  # p(z | tree, alpha, theta, sigma, sigmae, g0) x N(g0 | g0Prior$mean,
  # g0Prior$var). 
  loglik <- function(par, pruneInfo, memo = NULL) {
    if(parDigits >= 0) {
      par <- as.vector(round(par, parDigits))
    }
    
    atsseg0 <- spec$parMapping(par)
    
    if(useCpp & usempfr <= 0) {
      #if(verbose) 
      #  cat('.')
      # no use of high-precision floating point ops, so we call the faster C++ 
      # likelihood implementation
      val <- likPOUMM_lowLevelFun(
        pruneInfo$integrator, 
        alpha = atsseg0[1], theta = atsseg0[2],
        sigma = atsseg0[3], sigmae = atsseg0[4], g0 = atsseg0[5], 
        g0Prior = spec$g0Prior, log = TRUE)
    } else {
      atsseg0.list <- as.list(atsseg0)
      atsseg0[[4]] <- sqrt(atsseg0.list[[4]]^2+se^2)
      val <- 
        do.call(likPOUMMGivenTreeVTips, c(list(z, tree), atsseg0.list, list(
          g0Prior = spec$g0Prior, log = TRUE, pruneInfo = pruneInfo, 
          usempfr = usempfr, ...)))
    }
    
    if(is.na(val) | is.infinite(val)) {
      val <- -1e100
      attr(val, "g0") <- atsseg0[5]
      attr(val, "g0LogPrior") <- NA
    } 
    
    if(!is.null(memo)) {
      valMemo <- mget('val', memo, ifnotfound = list(-Inf))$val
      valPrev <- mget('valPrev', memo, ifnotfound = list(-Inf))$valPrev
      valDelta <- mget('valDelta', memo, ifnotfound = list(NA))$valDelta  
    } else {
      valMemo <- valPrev <- valDelta <- NA
    }
    
    if(usempfr == 0 & 
       (!is.na(valMemo) & valMemo < val & 
        ((valMemo + 1 < val) | 
         (valPrev + .01 * abs(valPrev) < val)))) {
  
        if(interactive() && verbose) {
          print(par)
          cat('Rmpfr check on: val =', val, '; valDelta =', val - valPrev)
        }
        
        atsseg0.list <- as.list(atsseg0)
        atsseg0[[4]] <- sqrt(atsseg0.list[[4]]^2+se^2)
      
        val2 <- 
          do.call(likPOUMMGivenTreeVTips, c(list(z, tree), atsseg0.list, list(
            g0Prior = spec$g0Prior, log = TRUE, pruneInfo = pruneInfo, 
            usempfr = 2, 
            ...)))
        
        if(is.na(val2) | is.infinite(val2)) {
          val2 <- -1e100
          attr(val2, "g0") <- atsseg0[5]
          attr(val2, "g0LogPrior") <- NA
        } 
        
        if(interactive() && verbose && abs(val2 - val) > abs(val) * 0.0001) {
          cat(' ====>  Changed difference - after Rmpfr-check valDelta =', val2 - valPrev, '.\n')
        } 
        val <- val2
    }
    
    val
  }
  
  result$loglik <- loglik
  
  if(interactive() && verbose) {
    print('Performing ML-fit...')
  }
  
  # default Rmpfr strategy 
  defaultRmpfr <- usempfr == 0
  
  # disable Rmpfr the first time we call maxLikPOUMMGivenTreeVTips
  if(defaultRmpfr) {
    usempfr <- -.5
  }
  fitML <- do.call(
    maxLikPOUMMGivenTreeVTips, 
    c(list(loglik = loglik, verbose = verbose, debug = debug, 
           pruneInfo = pruneInfo), spec))
  
  
  if(defaultRmpfr) {
    if(interactive() && verbose) {
      cat('Checking the max-loglik value with Rmpfr, current: val = ', 
          fitML$value, ", par=(", toString(fitML$par), ")")
    }
    usempfr = 2
    valLoglikRmpfr <- loglik(fitML$par, pruneInfo)
    if(interactive() && verbose) {
      cat(' ===> New: ', valLoglikRmpfr)
    }
    if(fitML$value - valLoglikRmpfr > 1E-6 * pruneInfo$N) {
      if(interactive()) {
        cat('Significant difference with Rmpfr-checked log-likelihood. Repeating the ML-fit with enabled Rmpfr checks for every improved log-likelihood point. You can disable this numerical stability test by setting usempfr = -1. Thanks for your patience.')
      }
      
      
      usempfr <-  0
      fitML <- do.call(
        maxLikPOUMMGivenTreeVTips, 
        c(list(loglik = loglik, verbose = verbose, debug = debug, 
               pruneInfo = pruneInfo), spec))
    } else {
      usempfr <- 0
      if(interactive() && verbose) {
        cat( " ===> OK.")
      }
    }
  } 
  
  
  if(interactive() && verbose) {
    cat("max loglik from ML: \n")
    print(fitML$value)
    cat("parameters at max loglik from ML: \n")
    print(fitML$par)
  }
  
  result[['fitML']] <- fitML
  
  # by default, no better likelihood found by mcmc (to be changed later if needed)
  result[['MCMCBetterLik']] <- 0 
  
  if(doMCMC & spec$nSamplesMCMC > 0) {
    if(interactive() && verbose) {
      print('Performing MCMC-fit...')
    }
    
    fitMCMC <- do.call(
      mcmcPOUMMGivenPriorTreeVTips, 
      c(list(loglik = loglik, fitML = fitML, verbose = verbose, debug = debug, 
             pruneInfo = pruneInfo), spec))
    
    if(interactive() && verbose) {
      cat("max loglik from MCMCs: \n")
      print(fitMCMC$valueMaxLoglik)
      cat("parameters at max loglik from MCMCs: \n")
      print(fitMCMC$parMaxLoglik)
    }
    
    # if the max likelihood from the MCMC is bigger than the max likelihood 
    # from the ML-fit, it is likely that the ML fit got stuck in a local 
    # optimum, so we correct it.
    if(fitML$value < fitMCMC$valueMaxLoglik) {
      parInitML <- as.vector(fitMCMC$parMaxLoglik)
      names(parInitML) <- names(spec$parLower)
      
      if(all(c(parInitML >= spec$parLower, parInitML <= spec$parUpper))) {
        if(interactive()) {
          cat("The MCMC-fit found a better likelihood than the ML-fit. Performing ML-fit starting from the MCMC optimum.")
        }
        spec[["parInitML"]] <- parInitML
        
        if(defaultRmpfr) {
          usempfr <- -.5
        }
        fitML2 <- do.call(
          maxLikPOUMMGivenTreeVTips, 
          c(list(loglik = loglik, verbose = verbose, debug = debug, 
                 pruneInfo = pruneInfo), spec))
        
        
        if(defaultRmpfr) {
          if(interactive() && verbose) {
            cat('Checking the max-loglik value with Rmpfr, current: val = ', fitML2$value)
          }
          usempfr = 2
          valLoglikRmpfr <- loglik(fitML2$par, pruneInfo)
          if(interactive() && verbose) {
            cat(' ===> New: ', valLoglikRmpfr)
          }
          if(fitML2$value - valLoglikRmpfr > 1E-6 * pruneInfo$N) {
            if(interactive()) {
              cat('Significant difference with Rmpfr-checked log-likelihood. Repeating the ML-fit with enabled Rmpfr checks for improved log-likelihood points. \nYou can disable this numerical stability test by setting usempfr = -1. Thanks for your patience.\n')
            }
            usempfr <-  0
            fitML2 <- do.call(
              maxLikPOUMMGivenTreeVTips, 
              c(list(loglik = loglik, verbose = verbose, debug = debug, 
                     pruneInfo = pruneInfo), spec))
          } else {
            usempfr <- 0
            if(interactive() && verbose) {
              cat( " ===> OK.")
            }
          }
        } 
        
        result[['fitML']] <- fitML2
        result[['MCMCBetterLik']] <- 1 # within search-region
      } else {
        message <- "The MCMC-fit found a better likelihood outside of the search-region of the ML-fit."
        if(interactive()) {
          cat(message)
        }
        result[['MCMCBetterLik']] <- 2 # outside search-region
      }
    }
    
    result[['fitMCMC']] <- fitMCMC
  }
  class(result) <- c('POUMM', class(result))
  result
}


#' Extract maximum likelihood and degrees of freedom from a fitted POUMM model
#' @param object An object of class POUMM.
#' @param ... not used; included for compliance with generic function logLik.
#' @export
logLik.POUMM <- function(object, ...) {
  if("POUMM" %in% class(object)) {
    lik <- object$fitML$value
    if(!is.null(object$dof)) {
      attr(lik, "df") <- object$dof
    } else {
      attr(lik, "df") <- length(coef(object))
    }
    attr(lik, "nobs") <- object$N
    class(lik) <- "logLik"
    lik
  } else {
    stop("logLik.POUMM called on non POUMM-object.")
  }
}

#' Extract maximum likelihood fitted parameters (coefficients) from a fitted 
#' POUMM  model.
#' 
#' @param object An object of class POUMM.
#' @param mapped Logical indicating whether the standard POUMM parameters should
#' also be extracted.
#' @param ... Not used; added for compatibility with generic function coef.
#' @details The parameters extracted are the ones used as input to the model's
#' parMapping function. 
#' 
#' 
#' @return A named vector with the fitted parameters of the model.
#' 
#' @importFrom stats coef
#' 
#' @export
coef.POUMM <- function(object, mapped = FALSE, ...) {
  if("POUMM" %in% class(object)) {
    pars <- object$fitML$par
    g0 <- attr(object$fitML$value, "g0")
    if(mapped) {
      pars <- object$spec$parMapping(pars)
    } else {
      if(is.null(names(pars))) {
        names(pars) <- names(object$spec$parLower)
      }
    }
    if("g0" %in% names(pars)) {
      pars["g0"] <- g0
    }
      
    pars
  } else {
    stop("coef.POUMM called on non POUMM-object.")
  }
}

#' Extract maximum likelihood expected genotypic values at the tips of a tree,
#' to which a POUMM model has been previously fitted
#' @param object An object of class POUMM.
#' @param vCov A logical indicating whether a list with the genotypic values and their variance covariance matrix should be returned or only a vector of the genotypic values (default is FALSE).
#' @param ... Not used; added for compatibility with generic function fitted.
#' @return If vCov == TRUE, a list with elements g - the genotypic values and
#' vCov - the variance-covariance matrix of these values for the specific tree,
#' observed values z and POUMM ML-fit. If vCov == FALSE, only the vector of 
#' genotypic values corresponding to the tip-labels in the tree is returned.
#' 
#' @importFrom stats fitted
#' @export
fitted.POUMM <- function(object, vCov=FALSE, ...) {
  if("POUMM" %in% class(object)) {
    g0 <- coef.POUMM(object, mapped=TRUE)['g0']
    if(is.nan(g0)) {
      if(interactive()) {
        cat("Genotypic values cannot be inferred for g0=NaN; Read documentaton for parMapping in ?specifyPOUMM and use a parMapping function that sets a finite value or NA for g0.")
      }
    }
    p <- coef(object, mapped = TRUE)
    gList <- gPOUMM(object$pruneInfo$z, object$pruneInfo$tree, g0,
                    p["alpha"], p["theta"], p["sigma"], p["sigmae"])
    g = as.vector(gList$mu.g.poumm)
    names(g) <- object$pruneInfo$tree$tip.label
    if(vCov) {
      list(g = g, vCov = gList$V.g.poumm)
    } else {
      g
    }
  } else {
    stop("fitted.POUMM called on non POUMM-object.")
  }
}

#' Number of tips in a phylogenetic tree, POUMM has been fit on.
#' @param object An object of class POUMM.
#' @param ... Not used; added for compatibility with generic function nobs.
#' @return The number of tips in the tree, POUMM has been called on
#' 
#' @export
nobs.POUMM <- function(object, ...) {
  if("POUMM" %in% class(object)) {
    object$N
  } else {
    stop("nobs.POUMM called on a non POUMM-object.")
  }  
}

#' Extract maximum likelihood environmental contributions (residuals) at the tips of a tree, to which a POUMM model has been fitted.
#' @param object An object of class POUMM.
#' @param ... Not used; added for compatibility with generic function residuals.
#' @return The vector of e-values (residuals) corresponding to the tip-labels in
#'  the tree.
#'  
#' @importFrom stats fitted
#' 
#' @export
residuals.POUMM <- function(object, ...) {
  if("POUMM" %in% class(object)) {
    e <- object$pruneInfo$z - fitted(object)
    names(e) <- object$pruneInfo$tree$tip.label
    e
  } else {
    stop("residuals.POUMM called on a non POUMM-object.")
  }
}

#' Plots of a POUMM-fit
#' @param x An object of class POUMM.
#' @param type A character indicating the type of plot(s) to be generated.
#'   Defaults to "MCMC", resulting in a trace and density plot for the selected
#'   statistics (see argument stat).
#' @param doPlot Logical indicating whether a plot should be printed on the 
#'   currently active graphics device or whether to return a list of ggplot 
#'   objects for further processing. Defaults to TRUE.
#' @param interactive Logical indicating whether the user should press a key 
#'   before generating a next plot (when needed to display two or more plots).
#'   Defaults to TRUE. Meaningless if 
#'   doPlot = FALSE.
#' @param stat A character vector with the names of statistics to be plotted.
#'   These should be names from the stats-list (see argument statFunctions).
#'   Defaults to c("alpha", "theta", "sigma", "sigmae", "H2tMean", "H2tInf").
#' @param chain A vector of integers indicating the chains to be plotted. 
#' @param startMCMC,endMCMC,thinMCMC Integers used to extract a sample from the 
#'   MCMC-chain; passed to summary().
#' @param statFunctions Named list of statistics functions; passed to summary().
#' @param doZoomIn (type MCMC only) A logical value indicating whether the 
#'   produced plots should have a limitation on the x-axis according to an 
#'   expression set in zoomInFilter (see below). Default value is FALSE.
#' @param zoomInFilter A character string which evaluates as logical value. If 
#'   doZoomIn is set to TRUE, this filter is applied to each point in each MCMC
#'   chain and the data-point is filtered out if it evaluates to FALSE. This 
#'   allows to zoomIn the x-axis of density plots but should be used with caution,
#'   since filtering out points from the MCMC-sample can affect the kernel densities.
#'   Unfortunately, filtering out values is currently the only way to affect the
#'   limits of individual facets in ggplot2. The default value is a complicated 
#'   expression involving the HPD from all MCMC chains (normally one chain from the
#'   prior and 2 chains from the posterior):
#'   zoomInFilter = paste0("(stat %in% c('H2e','H2tMean','H2tInf','H2tMax') & ", 
#'   "(value >= 0 & value <= 1) ) |",
#'   "( !stat %in% c('H2e','H2tMean','H2tInf','H2tMax') & ",
#'   "(value <= median(HPDUpper) + 4 * (median(HPDUpper) - median(HPDLower)) &",
#'   "value >= median(HPDLower) - 4 * (median(HPDUpper) - median(HPDLower))))").
#'  The identifiers in this expression can be any
#'   column names found in a summary of a POUMM object.
#' @param prettyNames A logical indicating if greek letters and sub/superscripts 
#' should be used for the names of columns in the posterior density pairs-plot.
#' @param showUnivarDensityOnDiag A logical indicating if univariate density 
#' plots should be displaied on the main diagonal in the bivariate posterior plot.
#' Defaults to FALSE, in which case the column names are displayed on the diagonal.
#' @param ... not used, needed for consistency with the generic plot-function.
#'
#' @return If doPlot==FALSE, a named list containing a member called data of
#'   class data.table and several members of class ggplot.
#' @export
plot.POUMM <- 
  function(x, type=c("MCMC"), 
           doPlot = TRUE, interactive = TRUE,
           stat=c("alpha", "theta", "sigma", "sigmae", "g0", "H2tMean"),
           chain=NULL,
           startMCMC = NA, endMCMC = NA, thinMCMC = 1000, 
           statFunctions = statistics(x),
           doZoomIn = FALSE,
           zoomInFilter = paste0("(stat %in% c('H2e','H2tMean','H2tInf','H2tMax') & ", 
                  "(value >= 0 & value <= 1) ) |",
                  "( !stat %in% c('H2e','H2tMean','H2tInf','H2tMax') & ",
                  "(value <= median(HPDUpper) + 4 * (median(HPDUpper) - median(HPDLower)) &",
                  "value >= median(HPDLower) - 4 * (median(HPDUpper) - median(HPDLower))))"),
           prettyNames = TRUE, 
           showUnivarDensityOnDiag = FALSE,
           ...) {
  
  if("POUMM" %in% class(x)) {
    summ <- summary(x, mode = "expert", 
                    startMCMC = startMCMC, endMCMC = endMCMC, 
                    thinMCMC = thinMCMC, stats = statFunctions)
    plot(summ, doPlot = doPlot, interactive = interactive, 
         stat = stat, chain = chain, 
         doZoomIn = doZoomIn, zoomInFilter = zoomInFilter, 
         prettyNames = prettyNames, 
         showUnivarDensityOnDiag = showUnivarDensityOnDiag)
  } else {
    stop("plot.POUMM called on a non POUMM-object.")
  }
}

#' A vectorized expected covariance function for a given tree and a fitted POUMM
#' model
#' 
#' @param object an S3 object of class POUMM
#' @param corr logical indicating if an expected correlation function 
#' should be returned
#'   For non-ultrametric trees, usually the mean root-tip distance is used.
#' 
#' @return a function of three numerical parameters:
#'  tau - the phylogenetic distance between a two tips;
#'  tanc - the distance from the root to their most recent common ancestor.
#'  t - the root-tip distance (assuming that the two tips are at equal distance from the root)
#' 
#' @export
covFunPOUMM <- function(object, corr=FALSE) {
  if("POUMM" %in% class(object)) {
    function(tau, tanc, t) {
             #t = mean(nodeTimes(object$pruneInfo$tree, 
            #                    tipsOnly = TRUE))) {
      par <- object$spec$parMapping(coef(object))
      covPOUMM(par['alpha'], par['sigma'], par['sigmae'], 
               t = t,
               tau = tau, 
               tanc = tanc, corr = corr)
    } 
  } else {
    stop("In covFunPOUMM: object should be of S3 class POUMM.")
  }
}

#' A vectorized function returning HPD intervals of the expected covariance for 
#' a given tree and a fitted POUMM model
#'
#' @param object an S3 object of class POUMM
#' @param prob a Numerical between 0 and 1
#' @param corr logical indicating if an expected correlation HPD interval 
#' function should be returned.
#' @param ... additional parameters passed to summary.POUMM
#' 
#' @return a function of a numerica matrix x with 3 columns corresponding to
#'  tau, tanc and t (see covFunPOUMM). The function reteruns a two-colun matrix
#'  with the lower and upper limit of the HPD for each row in the input matrix.
#'  
#' @import coda
#' @export
covHPDFunPOUMM <- function(object, prob = .95, corr = FALSE, ...) {
  # avoid check note "no visible binding":
  stat <- NULL
  
  if("POUMM" %in% class(object)) {
    smm <- summary(object, mode="long", ...)
    setkey(smm, stat)
    mcmc_asse <- smm[list(c('alpha', 'sigma', 'sigmae')), mcmc]
    asse <- cbind(mcmc_asse[[1]], mcmc_asse[[2]], mcmc_asse[[3]])
    tMean <- mean(nodeTimes(object$pruneInfo$tree, tipsOnly = TRUE))
    
    function(x) {
      t(apply(x, 1, function(xi) {
        covars <- apply(asse, 1, function(asse) {
          covPOUMM(asse[1], asse[2], asse[3], 
                   tau = xi[1], tanc = xi[2], t = xi[3], corr = corr)
        })  
        covars_mcmc <- mcmc(covars, start = start(mcmc_asse[[1]]), 
                            end = end(mcmc_asse[[1]]), 
                            thin = thin(mcmc_asse[[1]]))  
        
        as.vector(HPDinterval(covars_mcmc, prob = prob))
      }))
    } 
  } else {
    stop("In corrFunPOUMM: object should be of S3 class POUMM.")
  }
}

#' Simulate a trait on a tree under a ML fit of the POUMM model
#' 
#' @description Use the maximum likelihood parameters of the model to simulate
#' trait values on a phylogenetic tree. 
#' 
#' @details This function is a shortcut to calling
#' \code{\link{rVNodesGivenTreePOUMM}}, which will map the inferred parameters 
#' of the model back to the original POUMM parameters alpha, theta, sigma, sigmae,
#' and g0. 
#' 
#' @param object an S3 object of class POUMM
#' @param tree a phylo object. If NULL (default) the trait is simulated on the
#' tree, on which the POUMM object has been fit.
#' 
#' 
#' @return a numerical vector containing the simulated trait value for each tip 
#' in the tree.
#' 
#' @seealso \code{\link{rVNodesGivenTreePOUMM}}
#' @export
simulateTrait <- function(object, tree = NULL) {
  if(is.null(tree)) {
    tree <- object$pruneInfo$tree
  }
  
  
  par <- object$spec$parMapping(coef(object))
  rVNodesGivenTreePOUMM(tree, z0 = par['g0'], alpha = par['alpha'], 
                        theta = par['theta'], sigma = par['sigma'],
                        par['sigmae'])
}