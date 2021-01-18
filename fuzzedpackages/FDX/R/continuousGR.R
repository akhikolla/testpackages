#'@name continuous.GR
#'@title Continuous Guo-Romano procedure
#'
#'@description
#'Apply the usual continuous [GR] procedure, with or without computing the
#'critical values, to a set of p-values. A non-adaptive version is available as
#'well.
#'
#'@details
#'\code{GR} and \code{NGR} are wrapper functions for \code{continuous.GR}. The
#'first one simply passes all its parameters to \code{continuous.GR} with
#'\code{adaptive = TRUE} and \code{NGR} does the same with
#'\code{adaptive = FALSE}.
#'
#'@seealso
#'\code{\link{kernel}}, \code{\link{FDX-package}}, \code{\link{continuous.LR}},
#'\code{\link{discrete.LR}}, \code{\link{discrete.GR}}, 
#'\code{\link{discrete.PB}}, \code{\link{weighted.LR}}, 
#'\code{\link{weighted.GR}}, \code{\link{weighted.PB}}
#'
#'@templateVar raw.pvalues TRUE
#'@templateVar pCDFlist FALSE
#'@templateVar alpha TRUE
#'@templateVar zeta TRUE
#'@templateVar direction FALSE
#'@templateVar adaptive TRUE
#'@templateVar critical.values TRUE
#'@templateVar exact FALSE
#'@templateVar pvalues FALSE
#'@templateVar sorted_pv FALSE
#'@templateVar stepUp FALSE
#'@templateVar support FALSE
#'@templateVar weights FALSE
#'@templateVar weighting.method FALSE
#'@template param 
#'
#'@template example
#'@examples
#'
#'GR.fast <- GR(raw.pvalues)
#'summary(GR.fast)
#'
#'GR.crit <- GR(raw.pvalues, critical.values = TRUE)
#'summary(GR.crit)
#'
#'NGR.fast <- NGR(raw.pvalues)
#'summary(NGR.fast)
#'
#'NGR.crit <- NGR(raw.pvalues, critical.values = TRUE)
#'summary(NGR.crit)
#'
#'@templateVar Critical.values TRUE
#'@templateVar Adaptive TRUE
#'@templateVar Weighting FALSE
#'@template return
#'
#'@importFrom stats qbeta pbinom
#'@importFrom DiscreteFDR match.pvals
#'@export
continuous.GR <- function(raw.pvalues, alpha = 0.05, zeta = 0.5, adaptive = TRUE, critical.values = FALSE){
  #--------------------------------------------
  #       check arguments
  #--------------------------------------------
  if(is.null(alpha) || is.na(alpha) || !is.numeric(alpha) || alpha < 0 || alpha > 1)
    stop("'alpha' must be a probability between 0 and 1!")
  
  if(is.null(zeta) || is.na(zeta) || !is.numeric(zeta) || zeta < 0 || zeta > 1)
    stop("'zeta' must be a probability between 0 and 1!")
  #--------------------------------------------
  #       number of tests and [alpha * m] + 1
  #--------------------------------------------
  m <- length(raw.pvalues)
  a <- floor(alpha * 1:m) + 1
  #--------------------------------------------
  #       Determine sort order and do sorting
  #--------------------------------------------
  o <- order(raw.pvalues)
  sorted.pvals <- raw.pvalues[o]
  #--------------------------------------------
  #        Compute [GR] or [NGR] significant p-values,
  #        their indices and the number of rejections
  #--------------------------------------------
  # compute transformed sorted p-values
  if(adaptive){
    y <- cummax(pbinom(a - 1, m - 1:m + a, pmin(1, sorted.pvals), lower.tail = FALSE)) # cummax(pbeta(pmin(1, sorted.pvals), a, m - 1:m + 1))
  }else{
    y <- cummax(pbinom(a - 1, m, pmin(1, sorted.pvals), lower.tail = FALSE)) # cummax(pbeta(pmin(1, sorted.pvals), a, m - a + 1))
  }
  # determine significant (transformed) p-values
  idx <- which(y > zeta)
  if(length(idx)){
    m.rej <- min(idx) - 1
    if (m.rej){
      # determine significant (observed) p-values in sorted.pvals
      idx <- which(raw.pvalues <= sorted.pvals[m.rej])
      pvec.rej <- raw.pvalues[idx]
    }else{
      idx <- numeric(0)
      pvec.rej <- numeric(0)
    }
  }else{
    m.rej <- m
    idx <- 1:m
    pvec.rej <- raw.pvalues
  }
  # find critical constants
  if(critical.values){
    if (adaptive){
      crit.constants <- qbeta(zeta, a, m - 1:m + 1)
    }else{
      crit.constants <- qbeta(zeta, a, m - a + 1)
    }
  }
  
  #--------------------------------------------
  #       Create output S3 object
  #--------------------------------------------
  output <- list(Rejected = pvec.rej, Indices = idx, Num.rejected = m.rej)
  
  # add adjusted p-values to output list
  ro <- order(o)
  output$Adjusted = y[ro]
  
  # add critical values to output list
  if(critical.values) output$Critical.values = crit.constants
  
  # include details of the used algorithm as strings
  alg <- "Continuous Guo-Romano procedure"
  output$Method <- if(!adaptive) paste("Non-Adaptive", alg) else alg
  output$FDP.threshold <- alpha
  output$Exceedance.probability <- zeta
  output$Adaptive <- adaptive
  
  # original test data
  output$Data <- list()
  output$Data$raw.pvalues <- raw.pvalues
  # object names of the data as strings
  output$Data$data.name <- deparse(substitute(raw.pvalues))
  
  class(output) <- "FDX"
  return(output)
}

#'@rdname continuous.GR
#'@export
GR <- function(raw.pvalues, alpha = 0.05, zeta = 0.5, critical.values = FALSE){
  return(continuous.GR(raw.pvalues, alpha, zeta, TRUE, critical.values))
}

#'@rdname continuous.GR
#'@export
NGR <- function(raw.pvalues, alpha = 0.05, zeta = 0.5, critical.values = FALSE){
  return(continuous.GR(raw.pvalues, alpha, zeta, FALSE, critical.values))
}