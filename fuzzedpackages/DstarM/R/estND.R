#' Estimate nondecision density
#'
#' @param res an object of class \code{D*M}.
#' @param tt optional timegrid if the nondecision density is to be estimated at a different grid than the model density.
#' @param data if \code{tt} is specified then the original dataset must be supplied too.
#' @param h Optional smoothing parameter to be used when estimating the nondecision model on a
#' different time grid than the decision model. If omitted, the smoothing parameter of the decision model
#' is used.
#' @param zp Zero padding the estimated nondecision density by this amount to avoid numerical artefacts.
#' @param upper.bound An upper bound for the nondecision density. Defaults to one.
#' Lowering this bound can increase estimation speed,
#' at the cost of assuming that the density of the nondecision distribution is zero past this value.
#' @param lower.bound A lower bound for the nondecision density. Defaults to zero.
#' Increasing this bound can increase estimation speed,
#' at the cost of assuming that the density of the nondecision distribution is zero past this value.

#' @param Optim a named list with identical arguments to \code{\link{DEoptim.control}}.
#' In addition, if \code{verbose} == TRUE \code{Optim$steptol} can be a vector, i.e. \code{c(200, 50, 10)} means:
#' Do 200 iterations then check for convergence, do 50 iterations then check for convergence,
#' check every 10 iterations for convergence until itermax is reached.
#' If there are multiple nondecision distributions to estimate, one can supply different estimation
#' parameters for every nondecision distribution by supplying Optim as a list of lists. Every sublists
#' then corresponds to parameters for one nondecision distribution and should consist of arguments for
#' \code{\link{DEoptim.control}}.
#' Defaults to \code{Optim = list(reltol = 1e-6, itermax = 1e4, steptol = 200, CR = .9, trace = 0)}.
#' @param verbose Numeric, should intermediate output be printed? Defaults to 1, higher values result in more progress output.
#' Estimation will speed up if set to 0. If nonzero, \code{Optim$trace} will be
#' forced to 0, hereby disabling the build in printing of \code{DEoptim}. To enable the
#' printing of \code{DEoptim}, set \code{verbose} to 0 and specify \code{Optim$trace}.
#' @param dist A matrix where columns represent nondecision distributions.
#' If this argument is supplied then the objective function will be evaluated in these values.
#' @param NDindex A vector containing indices of which nondecision distributions to estimate.
#' If omitted, all nondecision distributions that complement the results in \code{res} are estimated.
#' @param max A positive float which indicates the maximum height of the nondecision distribution.
#' If estimated nondecision distributions appear chopped of or have a lot of values at this \code{max}
#' value it is recommended to re-estimate the nondecision distributions with a higher max value. Increasing
#' the \code{max} value without reason will increase the size of the parameter space and slow the estimation
#' procedure.
#' @param useRcpp Logical, setting this to true will make use of an Rcpp implementation of the objective function.
#' This gains speed at the cost of flexibility.
#'
#' @details
#' When verbose is set to 1, the ETA is an estimated of the time it takes to execute ALL iterations.
#' Convergence can (and is usually) reached before then.
#'
#' @examples
#' # simulate data with three stimuli of different difficulty.
#' # this implies different drift rates across conditions.
#' # define a time grid. A more reasonable stepsize is .01; this is just for speed.
#' tt = seq(0, 5, .1)
#'pars = c(.8, 2, .5, .5, .5, # condition 1
#'         .8, 3, .5, .5, .5, # condition 2
#'         .8, 4, .5, .5, .5) # condition 3
#'pdfND = dbeta(tt, 10, 30)
#'# simulate data
#'dat = simData(n = 3e5, pars = pars, tt = tt, pdfND = pdfND)
#'# define restriction matrix
#'restr = matrix(1:5, 5, 3)
#'restr[2, 2:3] = 6:7 # allow drift rates to differ
#'# fix variance parameters
#'fixed = matrix(c('sz1', .5, 'sv1', .5), 2, 2)
#'\dontrun{
#'# Run D*M analysis
#'res = estDstarM(data = dat, tt = tt, restr = restr, fixed = fixed)
#'# Estimate nondecision density
#'resND = estND(res)
#'plot(resND)
#'lines(tt, pdfND, type = 'b', col = 2)
#'}



#' @export
# estimate nondecision distribution
estND <- function(res, tt = NULL, data = NULL, h = res$h, zp = 5, upper.bound = 1,
  lower.bound = 0, Optim = list(), verbose = TRUE, dist = NULL, NDindex,
  max = 100, useRcpp = TRUE) {
  # zp - zero padding: adding a number of zeros to avoid numerical
  # artefacts these will be forced to 0 after the estimation procedure
  stopifnot(is.DstarM.fitD(res))

  if (xor(is.null(tt), is.null(data))) {
    stop("Supply both a time grid and the data to calculate the nondecision model at a custom time grid. Only one was supplied.")
  }

  ncondition <- res$ncondition
  if (!any(is.null(tt), is.null(data))) {
    # necessary to recalculate time grids

    data <- getData(res[["formula"]], data, verbose = verbose)
    rtime <- data[["rtime"]]
    response <- data[["response"]]
    condition <- data[["condition"]]
    hasConditions <- data[["hasConditions"]]
    data <- data[["data"]]

    if (is.null(h)) {
      # for backwards compatibility
      h <- 1
    }
    by <- unique(zapsmall(diff(tt)))
    note <- errCheckData(data = data, tt = tt, h = h, by = by, rtime = rtime,
      response = response, condition = condition)

    ncondition <- length(unique(data$condition))  # get number of conditions

    if (ncondition != res$ncondition) {
      warning(sprintf("Number of conditions in supplied data (%d) does not match number of conditions in decision model (%d)",
        ncondition, res$ncondition), call. = FALSE, immediate. = TRUE)
    }
    # recalculate modelDist and g to the accuracy of tt
    stopifnot(length(unique(round(diff(tt), 10))) == 1)  # also copy to objective.wrapper
    mm <- matrix(0, ncondition * 2, ncondition)
    mm[1:dim(mm)[1L] + dim(mm)[1L] * rep(1:dim(mm)[2L] - 1, each = 2)] <- 1

    rt <- split(data[[rtime]], list(data[[response]], data[[condition]]))

    g <- getGhat(rt, tt, ncondition, mm = mm, by = by)
    # convolve data density with uniform kernel
    kernel <- rev(stats::dunif(tt, 0, 1))  # 1 moet h zijn
    for (i in 1:dim(g)[2L]) {
      g[, i] <- customConvolveO(as.double(g[, i]), kernel)[seq_along(tt)]
    }
    g[g < 0] <- 0
    g <- g %*% (diag(dim(g)[2L])/rep(apply(g %*% mm, 2, simpson, x = tt),
      each = 2))

    pars <- res$Bestvals[c(res$restr.mat)]  # extract all parameters
    dim(pars) <- dim(res$restr.mat)  # set proper dimensions
    pars.list <- unlist(apply(pars, 2, list), recursive = FALSE, use.names = FALSE)  # reshape to lists
    modelDist <- getPdf(pars.list, tt, DstarM = TRUE, oscPdf = FALSE,
      mm, fun.density = res$fun.density)  # get Voss pdf
  } else {
    # not necessary to recalculate time grids
    tt <- res$tt
    g <- res$g.hat
    modelDist <- res$modelDist
  }
  group <- groups(ncondition, res$splits)
  k <- rev(stats::dunif(tt, 0, 1))
  by <- round(tt[2] - tt[1], 10)

  if (missing(NDindex)) {
    NDindex <- 1:dim(group)[2L]
  } else {
    NDindex <- unique(NDindex)
    stopifnot(all(NDindex > 0), length(NDindex) <= length(1:dim(group)[2L]),
      NDindex %in% 1:dim(group)[2L])
  }
  # if there is 1 list and no sublists this is 0, otherwise the number of
  # sublists
  tmp <- Reduce("+", lapply(Optim, is.list))
  tmp <- ifelse(is.integer(tmp), tmp, 0)
  if (tmp == 0) {
    OptimAll <- errCheckOptim(Optim, values = c(1e-04, 10000, 200, 0.9,
      0, 0))
    OptimAll <- rep(list(OptimAll), length(NDindex))
  } else if (tmp == length(NDindex)) {
    OptimAll <- lapply(Optim, errCheckOptim, values = c(1e-04, 10000,
      200, 0.9, 0, 0))
  } else {
    stop("Optim must be eiter a list of lists with length(NDindex) or a list with arguments for DEoptim.control().",
      call. = FALSE)
  }
  if (!is.null(dist)) {
    dist <- as.matrix(dist)  # to allow vector input
    if (dim(dist)[2L] != dim(group)[2L]) {
      stop(sprintf("dist must be a matrix with at least %s columns. Each column corresponds to a nondecision distribution to evaluate the objective function in.",
        dim(group)[2L]))
    } else {
      objVals <- rep.int(NA, dim(group)[2L])
      names(objVals) <- 1:dim(group)[2L]
    }
  }
  if (is.vector(lower.bound) & is.vector(upper.bound)) {
    if (length(lower.bound) == 1) {
      lower.bound <- rep(lower.bound, length(NDindex))
    } else if (length(lower.bound) != length(NDindex)) {
      stop("lower.bound must have either length one or length(NDindex).")
    }
    if (length(upper.bound) == 1) {
      upper.bound <- rep(upper.bound, length(NDindex))
    } else if (length(upper.bound) != length(NDindex)) {
      stop("upper.bound must have either length one or length(NDindex).")
    }
  } else {
    stop("lower.bound and/ or upper.bound must a vector.")
  }
  res.r <- vector("list", length = length(NDindex))
  ttr <- vector("list", length = length(NDindex))
  for (j in NDindex) {
    # j may select 1, 4, 5
    jidx <- which(NDindex == j)  # jidx will be 1, 2, 3
    Optim <- OptimAll[[jidx]]
    nrep <- ceiling((Optim$itermax - sum(Optim$steptol[1:(length(Optim$steptol) -
      1)]))/Optim$steptol[length(Optim$steptol)])
    if (nrep <= 1) {
      nrep <- which.max(cumsum(Optim$steptol) * (cumsum(Optim$steptol) <=
        Optim$itermax))
    } else if (nrep > 1) {
      last <- length(Optim$steptol)
      Optim$steptol2 <- rep(0, length(Optim$steptol) + nrep)
      Optim$steptol2[1:(last)] <- Optim$steptol
      Optim$steptol2[(last + 1):length(Optim$steptol2)] <- Optim$steptol[last]
      Optim$steptol <- Optim$steptol2
      nrep <- length(Optim$steptol)
    }

    argsList <- list(tt = tt)

    ttr[[jidx]] <- seq.int(lower.bound[jidx], upper.bound[jidx], by)
    ttr2 <- seq.int(lower.bound[jidx], upper.bound[jidx] + by * zp, by)

    nPre <- length(seq(min(tt), min(ttr2), by)) - 1L
    nPost <- length(seq(max(ttr2), max(tt), by)) - 1L

    # select objective function
    if (useRcpp) {
      if (nPre == 0 && nPost == 0) {

        argsList$fn <- rObjC0

      } else if (nPre == 0) {

        argsList$fn <- rObjC1
        argsList$lenPost <- rep(0, nPost)

      } else if (nPost == 0) {

        argsList$fn <- rObjC2
        argsList$lenPre <- rep(0, nPre)

      } else {

        argsList$fn <- rObjC3
        argsList$lenPre <- rep(0, nPre)
        argsList$lenPost <- rep(0, nPost)
      }
    } else {

      argsList$fn <- r.obj
      argsList$lenPre <- rep(0, nPre)
      argsList$lenPost <- rep(0, nPost)

    }

    argsList$control <- DEoptim::DEoptim.control(itermax = Optim$itermax)
    sharedNames <- names(argsList$control)[(names(argsList$control) %in%
      names(Optim))]
    argsList$control[sharedNames] <- Optim[sharedNames]

    if (Optim$parallelType == 1) {
      argsList$control$packages <- c("DstarM", Optim$packages)
      argsList$control$parVar <- Optim$parVar
    } else if (Optim$parallelType == 2) {
      argsList$control$foreachArgs <- c(list("simpson", "customConvolveO",
        "chisq", "customApprox"), Optim$foreachArgs)
    }

    argsList$lower <- rep(0, length(ttr2))
    argsList$upper <- argsList$lower + max  # + 1e30 is arbitrary and needs to be something better!
    if (is.null(Optim$NP[jidx]) | anyNA(Optim$NP[jidx])) {
      argsList$control$NP <- 10 * length(argsList$lower)
    } else {
      argsList$control$NP <- Optim$NP[jidx]
    }

    # argsList$control$initialpop = NULL # moet iets anders zijn!
    ind <- stats::na.omit(group[, j])  # index to avoid having to write na.omit 1e6 times
    m <- modelDist[, ind]  # get model dists
    temp1 <- 1/ncondition * rowSums(g[, ind])  # eq from article
    temp2 <- apply(m, 2, customConvolveO, y = k)[seq_along(temp1), ]
    temp2 <- rowSums(temp2)/ncondition  # 1/I sum[K*Mp(O)]
    # normalize temp1 and temp 2 -- redundant? NO
    temp1 <- temp1/simpson(tt, temp1)
    temp2 <- by * temp2/simpson(tt, temp2)
    if (!useRcpp) {
      # reverse here instead of in objective function
      temp2 <- rev(temp2)
    }
    if (!is.null(dist)) {
      objVals[j] <- r.obj(r = dist[, jidx], tt = tt, a = temp1, bb = temp2,
        lenPre = numeric(0), lenPost = numeric(0))  #, by = by)
    } else {
      # do differential evolution update arguments
      argsList$a <- temp1
      argsList$bb <- temp2
      if (verbose) {
        # do DEoptim for steptol iterations with custom printing and custom
        # convergence checks
        oldval <- Inf
        nfeval <- 0
        argsList$control$itermax <- Optim$steptol[1L]
        argsList$control$storepopfrom <- Optim$steptol[1L] - 1L
        before <- Sys.time()
        for (i in 1:nrep) {
          out <- do.call(DEoptim::DEoptim, argsList)
          argsList$control$initialpop <- out$member$pop  # update population value
          argsList$control$itermax <- Optim$steptol[i + 1L] + 1L  # update population value
          argsList$control$storepopfrom <- Optim$steptol[i + 1L]
          nfeval <- nfeval + out$optim$nfeval
          newval <- out$optim$bestval
          if (verbose == 1) {

            diff <- Sys.time() - before
            currentIter <- sum(Optim$steptol[1:i])
            totIter <- sum(Optim$steptol)
            # based on dplyr::progress_estimated
            avg <- diff/currentIter
            eta <- avg * (totIter - currentIter)

            cat(sprintf("Iteration: %d | Rel. improvement: %.3g | ETA: %.2f %s \n",
                        currentIter, (oldval - newval)/newval, eta, attr(eta, "unit")))
          } else if (verbose == 2) {
            replicate(120, cat("\b"))
            cat(sprintf("\rEstimating nondecision distribution %s out of %s \nTotal iterations done: %s \nImprovement over last %s iterations: %10g \nObjective function value: %10g",
                        jidx, length(NDindex), sum(Optim$steptol[1:i]), Optim$steptol[i],
                        oldval - newval, newval))
          }
          # check for convergence
          if (oldval - newval < Optim$reltol * newval) {
            break
          } else {
            oldval <- newval
          }
        }
        out$nfeval <- nfeval
      } else {
        # do DEoptim with build in convergence
        argsList$control$itermax <- Optim$itermax
        argsList$control$steptol <- Optim$steptol[1]
        argsList$control$reltol <- Optim$reltol
        out <- do.call(DEoptim::DEoptim, argsList)
      }
      res.r[[jidx]] <- out
    }
  }
  if (is.null(res[["conditionNames"]])) {
    # to remove: backwards compatibility
    names(res.r) <- paste0("ND", 1:length(res.r))
  } else {
    u <- unique(res$splits)
    nms <- sapply(u, function(x, all, nms) paste(nms[which(all == x)],
      collapse = "_"), all = res$splits, nms = res$conditionNames)
    names(res.r) <- nms
  }
  if (is.null(dist)) {
    r.hat <- matrix(data = 0, nrow = length(tt), ncol = length(res.r))
    # store estimated pars in common elements of ttr and ttt (rounding
    # because floats) also normalize estimated density so the area is 1.
    for (i in 1:length(res.r)) {
      r.hat[round(tt, 10) %in% round(ttr[[i]], 10), i] <- res.r[[i]]$optim$bestmem[seq_along(ttr[[i]])]  #res.r[[i]]$par[seq_along(ttr)]
      r.hat[, i] <- r.hat[, i]/simpson(tt, r.hat[, i])
    }
    GlobalOptimizer <- res.r
  } else {
    GlobalOptimizer <- objVals[!is.na(objVals)]
    r.hat <- dist
  }

  # Calculate some descriptives for the estimated densities
  descriptives <- matrix(NA, nrow = 6, ncol = dim(r.hat)[2L], dimnames = list(c("25%",
    "50%", "75%", "mean", "var", "mode"), names(res.r)))
  tr <- try(estQdf(p = seq(0, 1, 0.25), x = tt, cdf = estCdf(r.hat))[2:4,
    ], silent = TRUE)
  if (!is.numeric(tr)) {
    warning("Failed to calculate basic quantiles for nondecision distribution. Estimating with a more narrow time grid may solve this problem.")
  } else {
    descriptives[1:3, ] <- tr  # 25% quant, median, 75-quantile
  }
  descriptives[4, ] <- apply(r.hat, 2, nth.momentS, x = tt)  # mean
  descriptives[5, ] <- apply(r.hat, 2, nth.cmomentS, x = tt, nth = 2)  # variance
  descriptives[6, ] <- tt[apply(r.hat, 2, which.max)]  # mode
  # collect output

  out <- list(r.hat = r.hat, tt = tt, ttr = ttr, GlobalOptimizer = GlobalOptimizer,
    zp = zp, descriptives = descriptives)
  class(out) <- "DstarM.fitND"
  if (verbose & is.null(dist)) {
    cat("\nAnalyses complete!")
  }
  return(out)
}

# objective function for nondecision retrieval * by
r.obj <- function(r, tt, a, bb, lenPre, lenPost) {
  bb2 <- customConvolveO(c(lenPre, r, lenPost), bb)[seq_along(a)]
  return(chisq(tt, a, bb2))
}


