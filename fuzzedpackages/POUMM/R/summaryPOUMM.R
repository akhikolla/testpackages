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


#' Summarize the results of a POUMM-fit
#' 
#' @param object a POUMM object returned by POUMM-function (see ?POUMM).
#' @param ... Not used, but declared for consistency with the generic method summary.
#' @param startMCMC,endMCMC integers indicating the range of the MCMC chains
#' to be used for the analysis (excluding the initial warm-up phase)
#' @param thinMCMC thinning interval of the MCMC chain to avoid strong 
#' autocorrelation between sampled elements;
#' @param stats a named list of functions of the form function(par) { number },
#' which are called for each sample of each mcmc chain in object. Defaults to 
#' a call of statistics(object) returning a list of statistics functions relevant for 
#' the object. See also statistics.
#' @param mode a character indicating the desired format of the returned summary 
#' as follows:
#' 'short' - a data.table with the ML and MCMC estimates of heritability, 
#' model parameters, root-value and other statistics. 
#' 'long' - same information as in 'short' but including also the samples, which
#' can be convenient for 
#' 
#' @import data.table
#' @import coda
#' @importFrom stats AIC
#' 
#' @export
summary.POUMM <- function(object, ...,
                          startMCMC = NA, endMCMC = NA, thinMCMC = 1000, 
                          stats = statistics(object),
                          mode = c('short', 'long', 'expert')) {
  
  # declare global variables to avoid CRAN CHECK NOTES "no visible binding":
  N <- MLE <- samplePriorMCMC <- HPD <- HPD50 <- ESS <- HPDUpperFiltered <- 
    HPDLowerFiltered <- value <- HPDUpper <- HPDLower <- it <- PostMean <- mcs <- 
    ESS <- nChains <- chain <- G.R. <- stat <- Mean <-  NULL
  
  mode <- tolower(mode)
  
  tipTimes <- nodeTimes(object$pruneInfo$tree, tipsOnly = TRUE)
  tMax <- max(tipTimes)
  tMean <- mean(tipTimes)
  
  parLower <- matrix(object$spec$parLower, nrow = 1)
  parUpper <- matrix(object$spec$parUpper, nrow = 1)
  parML <- matrix(object$fitML$par, nrow = 1)
  
  anlist <- lapply(seq_along(stats), function(i) {
    data.table(stat = names(stats)[i], MLE = stats[[i]](parML))
  })
  
  anlist <- c(anlist, list(
    data.table(stat = "logpost", MLE = NA),
    data.table(stat = "loglik", MLE = object$fitML$value),
    data.table(stat = "AIC", MLE = AIC(object)),
    data.table(stat = "AICc", MLE = AIC(object) + 2*object$dof*(object$dof+1)/(object$N-object$dof-1))
  ))
  
  an.ML <- rbindlist(anlist)
  an.ML[, N:=object$N]
  setcolorder(an.ML, c('stat', 'N', 'MLE'))
  
  if(!is.null(object$fitMCMC)) {
    if(is.na(startMCMC)) {
      startMCMC <- object$spec$nSamplesMCMC / 10
    } 
    if(is.na(endMCMC)) {
      endMCMC <- object$spec$nSamplesMCMC
    } 
    
    anlist <- lapply(seq_along(stats), function(i) {
      analyseMCMCs(object$fitMCMC$chains, 
                   stat = stats[[i]], statName = names(stats)[i],
                   start = startMCMC, end = endMCMC, thinMCMC = thinMCMC, 
                   as.dt = TRUE)
    })
    
    anlist <- c(anlist, list(
      analyseMCMCs(object$fitMCMC$chains, 
                   stat=NULL, statName='logpost',
                   start = startMCMC, end=endMCMC, thinMCMC = thinMCMC, 
                   as.dt = TRUE),
      analyseMCMCs(object$fitMCMC$chains, 
                   stat = NULL, statName='loglik', 
                   start = startMCMC, end = endMCMC, thinMCMC = thinMCMC, 
                   as.dt = TRUE),
      analyseMCMCs(object$fitMCMC$chains, 
                   stat = NULL, statName='AIC', 
                   start = startMCMC, end = endMCMC, thinMCMC = thinMCMC, 
                   as.dt = TRUE, k = object$dof, N = object$N),
      analyseMCMCs(object$fitMCMC$chains, 
                   stat = NULL, statName='AICc', 
                   start = startMCMC, end = endMCMC, thinMCMC = thinMCMC, 
                   as.dt = TRUE, k = object$dof, N = object$N)
    ))
    
    an.MCMC <- rbindlist(anlist)
    an.MCMC[, samplePriorMCMC:=rep( c(object$spec$samplePriorMCMC, 
                                    rep(FALSE, object$spec$nChainsMCMC - 1)), length.out=.N) ]
    
    if(mode[1] != 'expert') {
      an.MCMC <- an.MCMC[, 
                         list(
                           PostMean = mean(unlist(Mean)),
                           HPD = list(colMeans(do.call(rbind, HPD))),
                           HPD50 = list(colMeans(do.call(rbind, HPD50))),
                           start = start(mcs), 
                           end = end(mcs), 
                           thin = thin(mcs),
                           ESS = sum(unlist(ESS)), 
                           G.R. = if(length(mcs)>1) { 
                             gelman.diag(mcs, autoburnin=FALSE)$psrf[1] 
                           } else {
                             as.double(NA)
                           },
                           nChains = length(mcs),
                           mcmc = mcmc.list(mcmc(do.call(rbind, mcs)))),
                         by=list(stat, samplePriorMCMC)]
    }
    
    if(mode[1] == 'short') {
      an.MCMC <- an.MCMC[
        samplePriorMCMC == FALSE, 
        list(stat, PostMean, HPD, ESS, G.R.)]
    } else if(mode[1] == 'long') {
      an.MCMC <- an.MCMC[
        samplePriorMCMC == FALSE, 
        list(stat, PostMean, HPD, HPD50, start, end, 
             thin = thinMCMC, ESS, G.R., nChains, mcmc)]
    } else if(mode[1] == 'expert') {
      an.MCMC <- an.MCMC[, list(stat, samplePriorMCMC, 
                                PostMean = Mean, HPD, HPD50, start, end, thin = thinMCMC,
                                ESS, mcmc = mcs, chain)]
    } else {
      warning(paste('mode should be one of "short", "long" or "expert", but was', mode[1], '.'))
    }
    
  } else {
    an.MCMC <- NULL
  }
  
  if(mode[1] %in% c('short', 'long')) {
    if(!is.null(an.ML) & !is.null(an.MCMC)) {
      res <- merge(an.ML, an.MCMC, by = 'stat', all = TRUE, sort = FALSE)
      res[sapply(HPD, is.null), HPD:=list(list(as.double(c(NA, NA))))]
      if(mode[1] == 'long') 
        res[sapply(HPD50, is.null), HPD50:=list(list(as.double(c(NA, NA))))]
    } else if(!is.null(an.ML)) {
      res <- an.ML
    } else if(!is.null(an.MCMC)) {
      res <- an.MCMC
    }
  } else {
    res <- list(spec = object$spec, startMCMC = startMCMC, endMCMC = endMCMC,
                thinMCMC = thinMCMC,
                ML = an.ML, MCMC = an.MCMC, 
                MCMCBetterLik = object$MCMCBetterLik)
  }
  
  class(res) <- c('summary.POUMM', class(res))
  res
}

generateStatisticFunG0 <- function(object) {
  function(par) {
    if( 'g0' %in% names(object$spec$parMapping(par)) ) {
      object$spec$parMapping(par)[, 'g0']
    } else {
      g0s <- try(
        apply(par, 1, function(par) {
          ll <- object$loglik(par, pruneInfo = object$pruneInfo)
          attr(ll, "g0")
        }), silent = TRUE)
      if(inherits(g0s, "try-error")) {
        rep(NA_real_, nrow(par))
      } else {
        g0s
      }
    }
  }
}


#' Plot a summary of a POUMM fit
#' @param x An object of class POUMM.
#' @param type A character indicating the type of plot(s) to be generated.
#'   Defaults to "MCMC", resulting in a trace and density plot for the selected
#'   statistics (see argument stat). Currently, only 'MCMC' type is supported.
#' @param doPlot Logical indicating whether a plot should be printed on the 
#'   currently active graphics device or whether only to return a list of plot- 
#'   objects for further processing. Defaults to TRUE.
#' @param stat A character vector with the names of statistics to be plotted.
#'   These should be names from the stats-list (see argument statFunctions).
#'   Defaults to c("alpha", "theta", "sigma", "sigmae", "H2tMean", "H2tInf").
#' @param chain A vector of integers indicating the chains to be plotted. 
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
#'    zoomInFilter = paste0("stat %in% c('H2e','H2tMean','H2tInf','H2tMax') |",
# "(value <= median(HPDUpper) + 4 * (median(HPDUpper) - median(HPDLower)) &",
#   "value >= median(HPDLower) - 4 * (median(HPDUpper) - median(HPDLower)))").
#'  The identifiers in this expression can be any
#'   column names found in a summary of a POUMM object.
#' @param palette A vector of colors (can be character strings) corresponding to the 
#' different chains (in their order 1 (prior), 2, 3). Defaults to c("#999999", 
#' "#0072B2", "#CC79A7", "#E69F00", "#D55E00", "#56B4E9", "#009E73", "#F0E442"),
#' which is a color-blind friendly.
#' @param prettyNames A logical indicating if greek letters and sub/superscripts 
#' should be used for the names of columns in the posterior density pairs-plot.
#' @param ... Not used; included for compatibility with the generic function plot.
#' 
#' @return If doPlot==TRUE, the function returns nothing and produces output on 
#' the current graphics device as a side-effect. Otherwise, the function returns
#' a list of plot-objects: traceplot and densplot.
#'  
#' @import ggplot2
#' @importFrom stats cor complete.cases
#' @import methods
#'  
#' @examples 
#' \dontrun{
#' library(POUMM)
#' 
#' set.seed(1)
#' 
#' N <- 1000
#' 
#' # create a random non-ultrametric tree of N tips
#' tree <- ape::rtree(N)  
#' 
#' # Simulate the evolution of a trait along the tree
#' z <- rVNodesGivenTreePOUMM(
#'   tree, g0 = 8, alpha = 1, theta = 4, sigma = 1.2, sigmae = .8)
#' 
#' fit <- POUMM(z[1:N], tree, spec = list(nSamplesMCMC = 4e5))
#' 
#' # Summarize the results from the fit in a table:
#' summary(fit)
#' 
#' # Create plots for some of the inferred parameters/statistics:
#' pl <- plot(fit, stat = c("alpha", "theta", "sigma", "sigmae", "H2tMean"), 
#'            doZoomIn = TRUE, 
#'            zoomInFilter = paste("!(stat %in% c('alpha', 'sigma', 'sigmae')) |",
#'                                 "(value >= 0 & value <= 8)"),
#'            doPlot = FALSE)
#' 
#' pl$traceplot
#' pl$densplot
#' }
#' 
#' @export
plot.summary.POUMM <- function(
  x, type = c("MCMC"), 
  doPlot = TRUE, 
  stat = c("alpha", "theta", "sigma", "sigmae", "g0", "H2tMean"),
  chain = NULL,
  doZoomIn = FALSE,
  zoomInFilter = paste0("(stat %in% c('H2e','H2tMean','H2tInf','H2tMax') & ", 
                        "(value >= 0 & value <= 1) ) |",
                        "( !stat %in% c('H2e','H2tMean','H2tInf','H2tMax') & ",
                        "(value <= median(HPDUpper) + 4 * (median(HPDUpper) - median(HPDLower)) &",
                        "value >= median(HPDLower) - 4 * (median(HPDUpper) - median(HPDLower))))"), 
  palette = c("#999999", "#0072B2", "#CC79A7", "#E69F00", "#D55E00", "#56B4E9", "#009E73", "#F0E442"), 
  prettyNames = TRUE,
  ...) {
  
  # declare global variables to avoid CRAN CHECK NOTES "no visible binding":
  N <- MLE <- samplePriorMCMC <- HPD <- HPD50 <- ESS <- HPDUpperFiltered <- 
    HPDLowerFiltered <- value <- HPDUpper <- HPDLower <- it <- PostMean <- mcs <- 
    ESS <- nChains <- chain <- G.R. <- stat2 <- statFactor <- NULL
  
  if(inherits(x, "summary.POUMM") && !is.null(x$MCMC)) {
    .stat <- stat
    .chain <- chain
    
    data <- merge(x$ML, x$MCMC, by = "stat")
    
    data <- data[
      { 
        if(!is.null(.stat)) {stat %in% .stat} else rep(TRUE, .N)
      } & { 
        if(!is.null(.chain)) {chain %in% .chain} else rep(TRUE, .N)
      }]
    
    setkey(data, stat)
    
    data <- data[list(.stat)]
    
    data <- data[{ 
      if(!is.null(.stat)) {stat %in% .stat} else rep(TRUE, .N)
    } & {
      if(!is.null(.chain)) {chain %in% .chain} else rep(TRUE, .N)
    }, list(
      N, 
      MLE,
      samplePriorMCMC,
      HPDLower = sapply(HPD, function(.) .[1]),
      HPDUpper = sapply(HPD, function(.) .[2]),
      HPD50Lower = sapply(HPD50, function(.) .[1]),
      HPD50Upper = sapply(HPD50, function(.) .[2]), 
      ESS,
      value = unlist(mcmc),
      it = seq(x$startMCMC, by = x$thinMCMC, along.with = mcmc[[1]])),

    by = list(stat = factor(stat), chain = factor(chain))]
    
    if(doZoomIn) {
      data[, stat2:=stat]
      
      data <- data[, {
        .SD[eval(parse(text = zoomInFilter))]
        }, by = stat2]
      
      data[, stat2:=NULL]
    }
    
    data[, HPDUpperFiltered:=min(max(value), unique(HPDUpper)), 
         list(stat = factor(stat), chain = factor(chain))]
    
    data[, HPDLowerFiltered:=max(min(value), unique(HPDLower)), 
         list(stat = factor(stat), chain = factor(chain))]
    
    .availStats <- data[, as.character(unique(stat))]
    
    statFactorLabels <- if(prettyNames) {
      prettifyNames(.availStats) 
    } else {
      .availStats
    }
    data[, statFactor:=factor(stat, levels = .availStats, labels = statFactorLabels)]
    
    .stat <- .availStats[1]
    dtm <- data[stat == .stat, list(stat, it, chain, value)]
    dtm[, (.stat) := value]
    dtm[, c("stat", "value") := NULL]
    
    for(.stat in .availStats[-1]) {
      dtm2 <- 
        data[stat == .stat, 
             eval(parse(text=paste0("list(it, chain, ", .stat, ' = value)')))]
      dtm <- merge(dtm, dtm2, by=c("it", "chain"), all=TRUE)
    }
    
    
    names(palette) <- as.character(seq_along(palette))
    
    my_ggplot <- function(...) ggplot(...) + 
      scale_color_manual(values = palette) +
      scale_fill_manual(values = palette)
    
    if(type == "MCMC") {
      traceplot <- my_ggplot(data) + 
        geom_line(aes(x=it, y=value, col = chain)) + 
        facet_wrap(~statFactor, 
                   scales = "free", 
                   labeller = if(prettyNames) "label_parsed" else "label_value")
      
      
      densplot <- my_ggplot(data) + 
        geom_density(aes(x=value, fill = chain, col = chain), alpha=0.5) + 
        geom_segment(aes(x=HPDLowerFiltered, xend=HPDUpperFiltered, 
                         y=0, yend=0, col = chain)) +
        geom_point(aes(x=MLE, y=0)) +
        facet_wrap(~statFactor, 
                   scales = "free", labeller= if(prettyNames) "label_parsed" else "label_value")
      
      if(doPlot) {
        print(traceplot) 
        if(interactive()) {
          print("Press Enter to see a univariate posterior density plot")
          scan("", what = "character", nlines = 1)
        }
        print(densplot) 
      } else {
        list(traceplot = traceplot, densplot = densplot)  
      }
    }
  } else {
    stop("plot.summary.POUMM called on a non summary.POUMM-object or a missing MCMC element. Verify that summary.POUMM has been called with mode = 'expert'")
  }
}

prettifyNames <- function(names) {
  prettyNames <- c(alpha = "alpha",
                   theta = "theta",
                   g0 = "g[0]", 
                   sigma = "sigma",
                   sigmae = "sigma[e]",
                   H2tMean = "H[bar(t)]^2",
                   H2e = "H[e]^2",
                   H2tInf = "H[infinity]^2")
  sapply(names, function(n) {pn <- prettyNames[n]; if(!is.na(pn)) pn else n}, USE.NAMES = FALSE)
}

#' Extract statistics from sampled or inferred parameters of a 
#' POUMM fit
#' @param object An object of class "POUMM".
#'
#' @details This is a generic method.
#' @export
statistics <- function(object) {
  UseMethod('statistics')
}

#' @describeIn statistics Relevant statistics from the sampled parameters of a
#'   POUMM fit
#'  
#' @export
statistics.POUMM <- function(object) {
  listPar <- sapply(seq_along(object$spec$parLower), function(i) {
    name <- names(object$spec$parLower)[i]
    stat <- eval(
      parse(text=paste0("list(", name, " = function(par) par[, ", i , "])"))
      )
  })
  listOtherStats <- list(
    H2e = function(par) H2e(z = object$pruneInfo$z,
                            sigmae = object$spec$parMapping(par)[, 'sigmae']),
    
    H2tInf = function(par) H2(alpha = object$spec$parMapping(par)[, 'alpha'],
                              sigma = object$spec$parMapping(par)[, 'sigma'],
                              sigmae = object$spec$parMapping(par)[, 'sigmae'],
                              t = Inf),
    
    H2tMax = function(par) H2(alpha = object$spec$parMapping(par)[, 'alpha'],
                              sigma = object$spec$parMapping(par)[, 'sigma'],
                              sigmae = object$spec$parMapping(par)[, 'sigmae'],
                              t = object$tMax),
    
    H2tMean = function(par) H2(alpha = object$spec$parMapping(par)[, 'alpha'],
                               sigma = object$spec$parMapping(par)[, 'sigma'],
                               sigmae = object$spec$parMapping(par)[, 'sigmae'],
                               t = object$tMean),
    
    alpha = function(par) object$spec$parMapping(par)[, 'alpha'],
    theta = function(par) object$spec$parMapping(par)[, 'theta'],
    sigma = function(par) object$spec$parMapping(par)[, 'sigma'],
    sigmae = function(par) object$spec$parMapping(par)[, 'sigmae'],
    g0 = generateStatisticFunG0(object),
    
    sigmaG2tMean = function(par) varOU(alpha = object$spec$parMapping(par)[, 'alpha'],
                                       sigma = object$spec$parMapping(par)[, 'sigma'],
                                       t = object$tMean),
    
    sigmaG2tMax = function(par) varOU(alpha = object$spec$parMapping(par)[, 'alpha'],
                                      sigma = object$spec$parMapping(par)[, 'sigma'],
                                      t = object$tMax),
    
    sigmaG2tInf = function(par) varOU(alpha = object$spec$parMapping(par)[, 'alpha'],
                                      sigma = object$spec$parMapping(par)[, 'sigma'],
                                      t = Inf))
  c(listPar, listOtherStats[setdiff(names(listOtherStats), names(listPar))])
}

