
mlegc <- function(y, x = NULL, locs, marginal, corr, effort = 1, longlat = FALSE, distscale = 1, method = 'GHK',
                    corrpar0 = NULL, ghkoptions = list(nrep = c(100,1000), reorder = FALSE, seed = 12345))
{
  cl <- match.call()
  if (!method %in% c("GQT", "GHK"))
    stop("'method' must be GQT or GHK.")

  if(!inherits(marginal, "marginal.gc"))  stop("'marginal' must be a function of the class marginal.gc")
  if(!inherits(corr, "corr.gc")) stop("'corr' must be a function of the class corr.gc")
  if(!isTRUE(marginal$discrete)) stop("marginals must be discrete")

  if (requireNamespace("numDeriv", quietly = TRUE)){
    if(method == 'GHK'){

      if(missing(ghkoptions)) ghkoptions <- list(nrep = c(100, 1000), reorder = FALSE, seed = 12345)
      if(is.null(ghkoptions[["nrep"]])) ghkoptions$nrep = c(100,1000)
      if(is.null(ghkoptions[["reorder"]])) ghkoptions$reorder = FALSE
      if(is.null(ghkoptions[["seed"]])) ghkoptions$seed = 12345
      y <- c(y)
      answer <- mleGHK(y = y, x = x, locs = locs, marginal = marginal, corr = corr, effort = effort,
                       longlat = longlat, distscale = distscale, corrpar0 = corrpar0, ghkoptions = ghkoptions)

      hessian <- numDeriv::hessian(func = likGHKXX, y = y, x = answer$MLE, XX = x, locs = locs, marginal = marginal, corr = corr, effort = effort,
                                   longlat = longlat, distscale = distscale, nrep = ghkoptions$nrep[length(ghkoptions$nrep)], seed = ghkoptions$seed)

    }else if(method == 'GQT'){
      answer <- mleGQT(y = y, x = x, locs = locs, marginal = marginal, corr = corr,
                       effort = effort, longlat = longlat, distscale = distscale, corrpar0 = corrpar0)

      hessian <- numDeriv::hessian(func = likGQTXX, y = y, x = answer$MLE, XX = x, locs = locs, marginal = marginal, corr = corr,
                                   effort = effort, longlat = longlat, distscale = distscale)
    }
  }else{
    stop("Please install {numDeriv} first!")
  }

  MIN <- .Machine$double.eps^0.25

  if(marginal$nod == 1 & answer$MLE[answer$kmarg] < MIN){
    warning("Estimated Standard Deviation of the Overdispersion Parameter may NOT Reliable Since
            its MLE is closed to the Boundary of Parameter Space!")
  }

  if(corr$nug == 1 & (answer$MLE[answer$par.df] < MIN | answer$MLE[answer$par.df] > 1-MIN)) {
    warning("Estimated Standard Deviation of the Nugget Parameter may NOT Reliable Since
            its MLE is closed to the Boundary of Parameter Space!")
  }

  if(corr$nug == 1 & answer$MLE[answer$par.df-1] < MIN) {
    warning("Estimated Standard Deviation of the Range Parameter may NOT Reliable Since
            its MLE is closed to the Boundary of Parameter Space!")
  }

  if(corr$nug == 0 & answer$MLE[answer$par.df] < MIN) {
    warning("Estimated Standard Deviation of the Range Parameter may NOT Reliable Since
            its MLE is closed to the Boundary of Parameter Space!")
  }

  answer[["hessian"]] <- hessian
  answer[["call"]] <- cl
  class(answer) <- c( "mlegc")
  return(answer)
  }









predgc <- function(obs.y, obs.x = NULL, obs.locs, pred.x = NULL, pred.locs, longlat = FALSE,
                     distscale = 1, marginal, corr, obs.effort = 1, pred.effort = 1, method = 'GHK',
                     estpar = NULL, corrpar0 = NULL, pred.interval = NULL, parallel = FALSE,
                     ghkoptions = list(nrep = c(100,1000), reorder = FALSE, seed = 12345),
                     paralleloptions = list(n.cores = 2, cluster.type = "SOCK"))
{
  if (!method %in% c("GQT", "GHK"))
    stop("'method' must be GQT or GHK.")
  if(!isTRUE(marginal$discrete)) stop("marginals must be discrete")

  if(is.null(pred.x)){
    nreg = 1
    }else{
    nreg <- ncol(pred.x)+1
    }
  nod <- marginal$nod
  ncorr <- corr$npar.cor
  Nestpar <- nreg + nod + ncorr

  if(!is.null(estpar)){
    if(length(estpar)!= Nestpar)
     stop("Number of the given parameter estimates does not match the model requirement, please check!")
  }

  if(parallel == FALSE & method == 'GHK'){
    answer <- predGHK(obs.y = obs.y, obs.x = obs.x, obs.locs = obs.locs, pred.x = pred.x, pred.locs = pred.locs,
                      longlat = longlat, distscale = distscale, marginal = marginal, corr = corr, obs.effort = obs.effort,
                      pred.effort = pred.effort, estpar = estpar, corrpar0 = corrpar0, pred.interval = pred.interval,
                      ghkoptions = ghkoptions)
  }else if(parallel == FALSE & method == 'GQT'){
    answer <- predGQT(obs.y = obs.y, obs.x = obs.x, obs.locs = obs.locs, pred.x = pred.x, pred.locs = pred.locs,
                      longlat = longlat, distscale = distscale, marginal = marginal, corr = corr, obs.effort = obs.effort,
                      pred.effort = pred.effort, estpar = estpar, corrpar0 = corrpar0, pred.interval = pred.interval)
  }else{

    #  if(0 != length(grep(paste("^package:", "snowfall", "$", sep=""), search())) == FALSE){
    #    stop("Please install {snowfall} first before using this function!")
    #}


    if(missing(paralleloptions)) paralleloptions <- list(n.cores = 2, cluster.type = "SOCK")
    if(is.null(paralleloptions[["n.cores"]])) paralleloptions$n.cores = 2
    if(is.null(paralleloptions[["cluster.type"]])) paralleloptions$cluster.type = "SOCK"

    if(method == 'GHK'){
      answer <- try(predGHK.sf(obs.y = obs.y, obs.x = obs.x, obs.locs = obs.locs, pred.x = pred.x, pred.locs = pred.locs,
                           longlat = longlat, distscale = distscale, marginal = marginal, corr = corr, obs.effort = obs.effort,
                           pred.effort = pred.effort, estpar = estpar,
                           corrpar0 = corrpar0, pred.interval = pred.interval,
                           n.cores = paralleloptions$n.cores,
                           cluster.type = paralleloptions$cluster.type,
                           ghkoptions = ghkoptions), silent = TRUE)
    }else{
      answer <- try(predGQT.sf(obs.y = obs.y, obs.x = obs.x, obs.locs = obs.locs, pred.x = pred.x, pred.locs = pred.locs,
                           longlat = longlat, distscale = distscale, marginal = marginal, corr = corr, obs.effort = obs.effort,
                           pred.effort = pred.effort, corrpar0 = corrpar0, estpar = estpar,
                           pred.interval = pred.interval,
                           n.cores = paralleloptions$n.cores, cluster.type = paralleloptions$cluster.type), silent = TRUE)
    }
    if(inherits(answer, "try-error")){
      stop("Please install and load the package {snowfall} first!")
    }
  }
  class(answer) <- c("predgc")
  return(answer)
}


