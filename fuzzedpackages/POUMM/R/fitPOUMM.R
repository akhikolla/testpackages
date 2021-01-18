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

# A utility function used to save the maximum point of f (used as wrapper for
# likelihod-functions).
memoiseMax <- function(f, par, memo, verbose, ...) {
  countMemo <- mget('count', envir = memo, ifnotfound = list(0))$count
  valMemo <- mget('val', envir = memo, ifnotfound = list(-Inf))$val
  valPrev <- mget('valPrev', envir = memo, ifnotfound = list(-Inf))$valPrev
  valDelta <- mget('valDelta', envir = memo, ifnotfound = list(NA))$valDelta
  parMemo <- mget('par', envir = memo, ifnotfound = list(NULL))$par
  
  assign("count", countMemo + 1, pos = memo)
  
  val <- f(par, memo = memo, ...)
  
  # if the returned val higher than the current max-value, store this val:
  if(valMemo < val) {
    assign('par', par, pos = memo)
    assign('val', val, pos = memo)
    assign('valDelta', val - valMemo, pos = memo)
    
    if(interactive() && verbose) {
      cat('\nCall ', countMemo, ': loglik on par=(', 
          toString(round(par, 6)), "): ", val, "\n", sep = "")
    }
  }
  
  # store the returned val and valDelta in memo
  assign('valPrev', val, pos = memo)
  
  val
}



#' Find a maximum likelihood fit of the POUMM model
#'
#' @param loglik function(par, memo, parFixedNoAlpha) 
#' @param pruneInfo a list-object returned by the pruneTree(tree, z) function.
#' @param parLower,parUpper Two named numeric vectors indicating the boundaries of 
#'   the search region. Default values are parLower = c(alpha = 0, theta = 0, 
#'   sigma = 0, sigmae = 0) and parUpper = c(alpha = 100, theta = 10, sigma = 20, 
#'   sigmae = 10).
#' @param parInitML A named vector (like parLower and parUpper) or a list of such
#'   vectors - starting points for optim.
#' @param control List of parameters passed on to optim, default 
#'   list(factr = 1e8, fnscale = -1), see ?optim. 
#' @param verbose A logical indicating whether to print informative messages on 
#'   the standard output.
#' @param debug A logical indicating whether to print debug messages 
#' (currently not implemented).
#' @param ... currently not used.
#'   
#' @return a list containing an element par and an element value as well as the 
#'   parameters passed
#'   
#' @importFrom stats optim runif
maxLikPOUMMGivenTreeVTips <- function(
  loglik, pruneInfo, 
  parLower, parUpper, parInitML = NULL, 
  control=list(factr = 1e8, fnscale = -1), 
  verbose = FALSE, debug = FALSE, ...) {
  if(is.null(parLower) | is.null(parUpper)) {
    stop("parLower and/or parUpper are NULL.")
  }
  
  
  if(is.null(parInitML)) {
    matParInitML <- rbind(
      parLower + 0.05 * (parUpper - parLower),
      parLower + 0.3 * (parUpper - parLower), 
      0.5 * (parLower + parUpper),
      parLower + 0.7 * (parUpper - parLower),
      parLower + 0.95 * (parUpper - parLower)
    )
    
    permuts_1to5 <- rep(
      c(4, 1, 3, 1, 2, 1, 2, 5, 4, 2, 5, 5, 2, 1, 3, 1, 3, 2, 3, 3, 4, 3, 5, 5, 
        5, 5, 4, 3, 2, 1, 4, 2, 3, 1, 1, 3, 4, 5, 2, 3, 1, 1, 4, 4, 2, 5, 4, 2, 
        4, 5), 2)
    
    
    matParInitML <- sapply(1:ncol(matParInitML), function(j) {
      matParInitML[, j][permuts_1to5[j:(j+50-1)]]
    })
      
    listParInitML <- lapply(1:nrow(matParInitML), function(i) matParInitML[i,])
    
  } else if(is.list(parInitML)) {
    listParInitML <- parInitML
  } else {
    listParInitML <- list(parInitML)
  }
  
  
  memoMaxLoglik <- new.env()
  fnForOptim <- function(par) {
    ll <- memoiseMax(
      loglik, par = par, pruneInfo = pruneInfo,  memo = memoMaxLoglik, verbose = verbose)
    g0LogPrior <- attr(ll, "g0LogPrior")
    
    if(is.finite(g0LogPrior)) {
      ll + g0LogPrior
    } else {
      ll
    }
  }
  
  for(iOptimTry in seq_along(listParInitML)) {
    parInitML <- listParInitML[[iOptimTry]]
    if(!(all(parInitML >= parLower) & all(parInitML <= parUpper))) {
      stop(paste0("All parameters in parInitML should be between \n parLower=", 
                  toString(round(parLower, 6)), " and \n parUpper=", 
                  toString(round(parUpper, 6)), ", but were \n ", 
                  toString(round(parInitML, 6)), "."))
    }
    
    if(is.null(control)) {
      control <- list()
    }
    
    # ensure that optim does maximization.
    control$fnscale <- -1
    
    if(length(parInitML) > 0) {
      if(interactive() && verbose) {
        cat("Call to optim no.", iOptimTry, ": starting from ", 
            toString(round(parInitML, 6)), "\n", 
            "parLower = ", toString(round(parLower, 6)), "\n",
            "parUpper = ", toString(round(parUpper, 6)), "\n")
      }
      res <- optim(fn = fnForOptim, 
                   par = parInitML, lower = parLower, upper = parUpper, 
                   method = 'L-BFGS-B', control = control)
    } else {
      stop("ML-fit: parInitML has zero length.")
    }
  }
  
  maxPar <- get("par", pos = memoMaxLoglik)
  maxValue <- get("val", pos = memoMaxLoglik)
  callCount <- get("count", pos = memoMaxLoglik)
  
  list(par = maxPar,
       value = maxValue, 
       g0 = attr(maxValue, which="g0", exact = TRUE), 
       g0LogPrior = attr(maxValue, which="g0LogPrior", exact = TRUE), 
       count = callCount, parLower=parLower, parUpper=parUpper, 
       control=control)
}


#' @importFrom stats window start end
#' @importFrom adaptMCMC convert.to.coda
convertToMCMC <- function(obj, thinMCMC=1) {
  
  codaobj <- convert.to.coda(obj)
  
  winmcmc <- window(codaobj, start = start(codaobj), end = end(codaobj), 
                    thin = thinMCMC)
  
  log.p <- obj$log.p
  
  log.p <- window(mcmc(log.p, 
                       start = start(codaobj), 
                       end = end(codaobj), 
                       thin = 1), 
                  start = start(codaobj), end = end(codaobj), thin = thinMCMC)

  list(
    mcmc = winmcmc, log.p = log.p, 
    accRateMCMC = obj$acceptance.rate, 
    adaption = obj$adaption, n.sample = obj$n.sample, cov.jump = obj$cov.jump,
    # This is causing the saved fit objects to be abnormally big.
    # Try not using it. 
    #sampling.parameters = obj$sampling.parameters, 
    scale.start = obj$scale.start)
}

#' Extract data from an MCMC chain
#' This is an internal function.
#' @param chains,stat,statName,start,end,thinMCMC,as.dt,k,N,... internal use.
#' 
#' @importFrom coda as.mcmc.list mcmc effectiveSize gelman.diag HPDinterval thin
#' 
analyseMCMCs <- function(chains, stat=NULL, statName="logpost", 
                         start, end, thinMCMC, as.dt=FALSE, k = NA, N = NA, ...) {
  mcmcs <- lapply(chains, function(obj) {
    if(is.null(obj$mcmc)) 
      convertToMCMC(obj, thinMCMC)
    else
      obj
  })
  
  log.ps <- window(as.mcmc.list(lapply(mcmcs, function(mc) { mc$log.p })),
                   start=start, end=end, thin=thinMCMC)
  
  if(statName == 'logpost') {
    mcs <- as.mcmc.list(
      lapply(seq_along(log.ps), function(i) {
          if(is.matrix(log.ps[[i]])) {
            log.ps[[i]][, 1, drop = FALSE]
          } else {
            mcmc(matrix(log.ps[[i]], length(log.ps[[i]])), 
                 start = start(log.ps[[i]]),
                 end = end(log.ps[[i]]), thin = thin(log.ps[[i]]))
                  
          }
        }))
    names(mcs) <- names(chains)
  } else if(statName == 'loglik') {
    mcs <- as.mcmc.list(
      lapply(seq_along(log.ps), function(i) {
        if(is.matrix(log.ps[[i]])) {
          log.ps[[i]][, 2, drop = FALSE]
        } else {
          mcmc(matrix(as.double(NA), length(log.ps[[i]])), 
               start = start(log.ps[[i]]),
               end = end(log.ps[[i]]), thin = thin(log.ps[[i]]))
        }
      }))
    names(mcs) <- names(chains)
  } else if(statName == 'AIC') {
    mcs <- as.mcmc.list(
      lapply(seq_along(log.ps), function(i) {
        if(is.matrix(log.ps[[i]])) {
          logl <- log.ps[[i]][, 2, drop = FALSE]
          2*k - 2*logl  
        } else {
          mcmc(matrix(as.double(NA), length(log.ps[[i]])), 
               start = start(log.ps[[i]]),
               end = end(log.ps[[i]]), thin = thin(log.ps[[i]]))
        }
      }))
    names(mcs) <- names(chains)
  } else if(statName == 'AICc') {
    mcs <- as.mcmc.list(
      lapply(seq_along(log.ps), function(i) {
        if(is.matrix(log.ps[[i]])) {
          logl <- log.ps[[i]][, 2, drop = FALSE]
          2*k - 2*logl + 2*k*(k+1)/(N-k-1)  
        } else {
          mcmc(matrix(as.double(NA), length(log.ps[[i]])), 
               start = start(log.ps[[i]]),
               end = end(log.ps[[i]]), thin = thin(log.ps[[i]]))
        }
        
      }))
    names(mcs) <- names(chains)
  # } else if(statName == 'g0') {
  #   mcs <- as.mcmc.list(
  #     lapply(seq_along(log.ps), function(i) {
  #       if(is.matrix(log.ps[[i]])) {
  #         log.ps[[i]][, 3, drop = FALSE]  
  #       } else {
  #         mcmc(matrix(as.double(NA), length(log.ps[[i]])), 
  #              start = start(log.ps[[i]]),
  #              end = end(log.ps[[i]]), thin = thin(log.ps[[i]]))
  #       }
  #     }))
  #   names(mcs) <- names(chains)
  # } 
    } else {
    mcs <- as.mcmc.list(lapply(mcmcs, function(mc) {
      winmcmc <- window(mc$mcmc, start = start, end = end, thin = thinMCMC)
      data <- matrix(stat(winmcmc, ...), nrow = nrow(winmcmc), byrow = TRUE)
      mcmc(data, start = start(winmcmc), end = end(winmcmc), 
                 thin = thin(winmcmc))
    }))
    names(mcs) <- names(chains)
  }
  
  mcs.mat <- as.matrix(mcs)
  
  if(length(chains) > 1) {
    gelman <- gelman.diag(mcs, autoburnin=FALSE)
    gelman <- ifelse(ncol(mcs.mat) > 1, gelman$mpsrf, gelman$psrf[2]) 
  } else {
    gelman <- NULL
  }
  ESS <- try(lapply(mcs, effectiveSize), silent=TRUE)
  if(inherits(ESS, 'try-error')) {
    warning(paste("For ", statName, 
                  "calling coda::effectiveSize resulted in an error:", 
                  toString(ESS)))
    ESS <- 0
  }
    
  names(ESS) <- NULL
  
  HPD <- lapply(mcs, function(.) {
    int <- try(HPDinterval(.), silent = TRUE)
    if(inherits(int, 'try-error')) {
      int <- matrix(as.double(NA), nrow = ncol(as.matrix(.)), ncol = 2)
      colnames(int) <- c("lower", "upper")
      int
    } else {
      int
    }
  })
  
  HPD50 <- lapply(mcs, function(.) {
    int <- try(HPDinterval(., 0.5), silent = TRUE)
    if(inherits(int, 'try-error')) {
      int <- matrix(as.double(NA), nrow = ncol(as.matrix(.)), ncol = 2)
      colnames(int) <- c("lower", "upper")
      int
    } else {
      int
    }
  })

  Mean <- sapply(seq_along(mcs), function(i) {
    colMeans(mcs[[i]])
  })
  
  if(as.dt) {
    data.table(
      chain=seq_along(mcs), stat=statName, start=start, end=end, thinMCMC=thinMCMC, 
      ESS=ESS, Mean=Mean, HPD=HPD, HPD50=HPD50, mcs=mcs)
  } else {
    l <- list(stat=statName, start=start, end=end, thinMCMC=thinMCMC, 
              ESS=ESS, Mean=Mean, HPD=HPD, HPD50=HPD50, G.R.=gelman, mcs=mcs, log.ps=log.ps)
    
    # without this we are running into trouble:
    names(l) <- c('stat', 'start', 'end', 'thinMCMC',
                  'ESS', Mean=Mean, 'HPD', 'HPD50', 'G.R.', 'mcs', 'log.ps')    
    l
  }
}

#' MCMC-sampling from a posterior distribution of a P(OU)MM model given tree, 
#' values at the tips and a prior distribution
#' 
#' @param loglik a log-likelihood function.
#' @param fitML an object returned by the maxLikPOUMMGivenTreeVTips
#' @param parMapping a function(numeric-vector) transforming a sampled vector on
#'   the scale of the parameters alpha, theta, sigma, sigmae and g0.
#' @param parInitMCMC a function(chainNumber) returning the starting point of 
#'   the MCMC as a vector.
#' @param parPriorMCMC a function(numeric-vector) returning the log-prior of the
#'   supplied vector
#' @param parScaleMCMC numeric matrix indicating the initial jump-distribution 
#'   matrix
#' @param nSamplesMCMC integer indicating how many iterations should the 
#'   mcmc-chain contain
#' @param nAdaptMCMC integer indicating how many initial iterations should be 
#'   used for adaptation of the jump-distribution matrix
#' @param thinMCMC integer indicating the thinning interval of the mcmc-chain
#' @param accRateMCMC (MCMC) numeric between 0 and 1 indicating the target 
#'   acceptance rate Passed on to adaptMCMC::MCMC.
#' @param gammaMCMC (MCMC) controls the speed of adaption. Should be in the interval (0.5,1]. A lower gammaMCMC leads to faster adaption. Passed on to 
#'   adaptMCMC::MCMC.
#' @param nChainsMCMC integer indicating the number of chains to run. Defaults 
#'   to 1.
#' @param samplePriorMCMC logical indicating if only the prior distribution 
#'   should be sampled. This can be useful to compare with mcmc-runs for an 
#'   overlap between prior and posterior distributions.
#' @param pruneInfo a list-object returned from the pruneTree(tree, z) function.
#' @param ... Additional arguments. Currently not used except for the following:
#'   If ... includes debug = TRUE, some debug messages will be written also
#'   outside of the call to loglik.
#' @param parallelMCMC Logical indicating if chains should be run in parallel.
#' @param verbose Logical indicating if some informal messages should be written
#'   during run. This parameter is passed to loglik.
#'   
#' @details Currently, this function calls the MCMC function from the adaptMCMC 
#'   package.
#'   
#' @return a list of coda objects
#' 
#' @importFrom foreach foreach %dopar% %do%
#' @importFrom adaptMCMC MCMC
#'   
#' @aliases mcmcPOUMMGivenPriorTreeVTips mcmc.poumm
mcmcPOUMMGivenPriorTreeVTips <- function(
  loglik, fitML = NULL,
  parMapping, parInitMCMC, parPriorMCMC, parScaleMCMC, 
  nSamplesMCMC, nAdaptMCMC, thinMCMC, accRateMCMC, 
  gammaMCMC, nChainsMCMC, samplePriorMCMC, pruneInfo, ..., 
  verbose = FALSE, parallelMCMC = FALSE) {
  
  samplePriorMCMC <- c(samplePriorMCMC, rep(FALSE, nChainsMCMC - 1))
  
  if(!is.null(list(...)$debug)) {
    debug <- list(...)$debug
  } else {
    debug <- FALSE
  }
  
  if(interactive() && debug) {
    cat("Using prior function:\n")
    print(parPriorMCMC) 
  }
  
  
  post <- function(par, memoMaxLoglik, chainNo, pruneInfo) {
    pr <- parPriorMCMC(par)
    if(debug) {
      cat("par:", toString(round(par, 6)), 'prior:', pr, ";")
    }
    if(is.nan(pr) | is.na(pr)) {
      warning(paste0("NA or NaN prior for ", toString(par)))
      -Inf
    } else if(is.infinite(pr)) {
      pr
    } else {
      if(samplePriorMCMC[chainNo]) {
        ll <- 0
        attr(ll, "g0") <- par[5]
        attr(ll, "g0LogPrior") <- NA
      } else {
        ll <- memoiseMax(
          loglik, par = par, pruneInfo = pruneInfo, 
          memo = memoMaxLoglik, verbose = verbose)
        
        if(interactive() && debug) {
          cat("par: ", toString(round(par, 6)), ";", "ll:", ll)
        }
      }
      if(interactive() && debug) {
        cat("\n")
      }
      
      pr + ll
    }
  }
  
  doMCMC <- function(i, pruneInfo) {
    memoMaxLoglik <- new.env()
    
    if(interactive() && verbose) {
      cat('Chain No:', i, ', nSamplesMCMC = ', nSamplesMCMC, ', initial state: \n')
    }
    
    init <- parInitMCMC(i, fitML)
    
    if(interactive() && verbose) {
      print(init)  
    }
    
    ch <- convertToMCMC(
      MCMC(
        post, n = nSamplesMCMC, init = init, scale = parScaleMCMC, adapt = nAdaptMCMC, 
        acc.rate = accRateMCMC, gamma = gammaMCMC, 
        memoMaxLoglik = memoMaxLoglik, chainNo = i, pruneInfo = pruneInfo), 
      thinMCMC = thinMCMC)
    
    if(interactive()) {
      print(paste("Finished chain no:", i))
    }
    
    ch$parScaleMCMC.start <- parScaleMCMC    
    
    ch$parMaxLoglik <- mget("par", envir = memoMaxLoglik, 
                            ifnotfound = list(NULL))$par
    ch$valueMaxLoglik <- mget("val", envir = memoMaxLoglik,
                              ifnotfound = list(-Inf))$val
    
    ch
  }
  
  if(parallelMCMC) {
    chains <- foreach(
      i = 1:nChainsMCMC, .packages = c('coda', 'POUMM', 'adaptMCMC')) %dopar% {
        # need to reinitialize the integrator since it is a C++ object
        # pruneInfo$integrator <- Integrator$new()
        # 
        # pruneInfo$integrator$setPruningInfo(
        #   pruneInfo$z,  pruneInfo$se, 
        #   pruneInfo$tree$edge,
        #   pruneInfo$tree$edge.length, 
        #   pruneInfo$M, pruneInfo$N, 
        #   pruneInfo$endingAt, 
        #   pruneInfo$nodesVector, pruneInfo$nodesIndex, 
        #   pruneInfo$unVector, pruneInfo$unIndex)
        # 
        pruneInfo$integrator <- POUMM_AbcPOUMM$new(pruneInfo$tree, pruneInfo$z[1:pruneInfo$N], 
                                             if(length(pruneInfo$se) != pruneInfo$N) {
                                               rep(pruneInfo$se[1], pruneInfo$N)
                                             } else {
                                               pruneInfo$se
                                             })
        
        chain <- try(doMCMC(i, pruneInfo = pruneInfo), silent = TRUE)
      }
  } else {
    chains <- foreach(
      i = 1:nChainsMCMC, .packages = c('coda', 'POUMM', 'adaptMCMC')) %do% {
        chain <- try(doMCMC(i, pruneInfo = pruneInfo), silent = TRUE)
      }
  }
  
  for(i in seq_along(chains)) {
    if(inherits(chains[[i]], "try-error")) {
      warning(paste0("Error in MCMC chain no ", i, ":", toString(chains[[i]])))
    }
  }
  # maximum log-likelihood from each chain. This is used to correct ML fit in
  # case it got stuck in a local optimum.
  maxLogliksFromChains <- sapply(chains, function(.) .$valueMaxLoglik)
  
  list(chains = chains, 
       post = post, 
       parMaxLoglik = chains[[which.max(maxLogliksFromChains)]]$parMaxLoglik,
       valueMaxLoglik = max(maxLogliksFromChains)
       )
}
