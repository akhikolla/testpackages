#'@name discrete.LR
#'
#'@title Discrete Lehmann-Romano procedure
#'
#'@description
#'Apply the [DLR] procedure, with or without computing the critical values, to
#'a set of p-values and their discrete support. Both step-down and step-up
#'procedures can be computed and non-adaptive versions are available as well.
#'
#'@details
#'\code{DLR} and \code{NDLR} are wrapper functions for \code{discrete.LR}.
#'The first one simply passes all its parameters to \code{discrete.LR} with
#'\code{adaptive = TRUE} and \code{NDLR} does the same with
#'\code{adaptive = FALSE}.
#'
#'@seealso
#'\code{\link{kernel}}, \code{\link{FDX-package}}, \code{\link{continuous.LR}},
#'\code{\link{continuous.GR}}, \code{\link{discrete.GR}}, 
#'\code{\link{discrete.PB}}, \code{\link{weighted.LR}},
#'\code{\link{weighted.GR}}, \code{\link{weighted.PB}}
#'
#'@templateVar raw.pvalues TRUE
#'@templateVar pCDFlist TRUE
#'@templateVar alpha TRUE
#'@templateVar zeta TRUE
#'@templateVar direction TRUE
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
#'@section References:
#' S. Döhler and E. Roquain (2019). Controlling False Discovery Exceedance for
#' Heterogeneous Tests.
#' \href{https://arxiv.org/abs/1912.04607v1}{arXiv:1912.04607v1}.
#'
#'@template example
#'@examples
#'
#'DLR.sd.fast <- DLR(raw.pvalues, pCDFlist)
#'summary(DLR.sd.fast)
#'DLR.su.fast <- DLR(raw.pvalues, pCDFlist, direction = "su")
#'summary(DLR.su.fast)
#'
#'DLR.sd.crit <- DLR(raw.pvalues, pCDFlist, critical.values = TRUE)
#'summary(DLR.sd.crit)
#'DLR.su.crit <- DLR(raw.pvalues, pCDFlist, direction = "su", critical.values = TRUE)
#'summary(DLR.su.crit)
#'
#'NDLR.sd.fast <- NDLR(raw.pvalues, pCDFlist)
#'summary(NDLR.sd.fast)
#'NDLR.su.fast <- NDLR(raw.pvalues, pCDFlist, direction = "su")
#'summary(NDLR.su.fast)
#'
#'NDLR.sd.crit <- NDLR(raw.pvalues, pCDFlist, critical.values = TRUE)
#'summary(NDLR.sd.crit)
#'NDLR.su.crit <- NDLR(raw.pvalues, pCDFlist, direction = "su", critical.values = TRUE)
#'summary(NDLR.su.crit)
#'
#'@templateVar Critical.values TRUE
#'@template return
#'
#'@importFrom DiscreteFDR match.pvals
#'@export
discrete.LR <- function(raw.pvalues, pCDFlist, alpha = 0.05, zeta = 0.5, direction = "sd", adaptive = TRUE, critical.values = FALSE){
  #--------------------------------------------
  #       check arguments
  #--------------------------------------------
  if(is.null(alpha) || is.na(alpha) || !is.numeric(alpha) || alpha < 0 || alpha > 1)
    stop("'alpha' must be a probability between 0 and 1!")
  
  if(is.null(zeta) || is.na(zeta) || !is.numeric(zeta) || zeta < 0 || zeta > 1)
    stop("'zeta' must be a probability between 0 and 1!")
  
  m <- as.integer(length(raw.pvalues))
  if(m != length(pCDFlist)) stop("The lengths of 'raw.pvalues' and 'pCDFlist' must be equal!") 
  #--------------------------------------------
  #       prepare p-values for processing
  #--------------------------------------------
  pvec <- match.pvals(pCDFlist, raw.pvalues)
  #--------------------------------------------
  #       Determine sort order and do sorting
  #--------------------------------------------
  o <- order(pvec)
  sorted.pvals <- pvec[o]
  #--------------------------------------------
  #       construct the vector of all values of all supports of the p-values
  #--------------------------------------------
  pv.list.all <- sort(unique(as.numeric(unlist(pCDFlist))))
  #--------------------------------------------
  #        Compute [HSU] or [HSD] significant p-values,
  #        their indices and the number of rejections
  #--------------------------------------------
  direction <- match.arg(direction, c("su", "sd"))
  if(direction == "su"){
    # SU case
    if(critical.values){
      # compute transformed support
      y <- kernel_DLR_crit(pCDFlist, pv.list.all, sorted.pvals, adaptive, alpha, zeta, TRUE)
      # find critical constants
      crit.constants <- y$crit.consts
      idx <- which(sorted.pvals <= crit.constants)
    }
    else{
      # compute transformed observed p-values
      y <- kernel_DLR_fast(pCDFlist, sorted.pvals, adaptive, alpha, TRUE, zeta, pv.list.all)
      # determine significant (transformed) p-values
      if(length(y)){
        idx <- which(y <= zeta * (floor(1:length(y) * alpha) + 1))
      }else{
        idx <- integer(0)
      }
    }
    if(length(idx)){
      m.rej <- max(idx)
      # determine significant (observed) p-values in sorted.pvals
      idx <- which(pvec <= sorted.pvals[m.rej]) 
      pvec.rej <- raw.pvalues[idx]
    }
    else{
      m.rej <- 0
      idx <- integer(0)
      pvec.rej <- numeric(0)
    }
  }
  else{
    # SD case
    if(critical.values){
      # compute transformed support
      y <- kernel_DLR_crit(pCDFlist, pv.list.all, sorted.pvals, adaptive, alpha, zeta, FALSE)
      # find critical constants
      crit.constants <- y$crit.consts
      idx <- which(sorted.pvals > crit.constants)
    }
    else{
      # compute transformed sorted p-values
      y <- kernel_DLR_fast(pCDFlist, sorted.pvals, adaptive, alpha, FALSE)
      # determine significant (transformed) p-values
      idx <- which(y > zeta * (floor(1:m * alpha) + 1)) 
    }
    if(length(idx)){
      m.rej <- min(idx) - 1
      if(m.rej){
        # determine significant (observed) p-values in sorted.pvals
        idx <- which(pvec <= sorted.pvals[m.rej])
        pvec.rej <- raw.pvalues[idx]
      }
      else{
        idx <- numeric(0)
        pvec.rej <- numeric(0)
      }
    }
    else{
      m.rej <- m
      idx <- 1:m
      pvec.rej <- raw.pvalues
    }
  }
  #--------------------------------------------
  #       Create output S3 object
  #--------------------------------------------
  output <- list(Rejected = pvec.rej, Indices = idx, Num.rejected = m.rej)
  if(direction == "sd"){
    if(critical.values){
      y <- y$pval.transf
    }
    # compute adjusted p-values
    pv.adj = cummax(pmin(y / (floor(alpha * 1:m) + 1), 1))
    # add adjusted p-values to output list
    ro <- order(o)
    output$Adjusted = pv.adj[ro]
  }
  # add critical values to output list
  if(critical.values) output$Critical.values = crit.constants
  
  # include details of the used algorithm as strings
  alg <- "Discrete Lehmann-Romano procedure"
  alg <- if(!adaptive) paste("Non-Adaptive", alg) else alg
  output$Method <- paste(alg, switch(direction, su = "(step-up)", sd = "(step-down)"))
  output$FDP.threshold <- alpha
  output$Exceedance.probability <- zeta
  output$Adaptive <- adaptive
  
  # original test data (often included, e.g. when using 'binom.test()')
  output$Data <- list()
  output$Data$raw.pvalues <- raw.pvalues
  output$Data$pCDFlist <- pCDFlist
  # object names of the data as strings
  output$Data$data.name <- paste(deparse(substitute(raw.pvalues)), "and", deparse(substitute(pCDFlist)))
  
  class(output) <- "FDX"
  return(output)
}

#'@rdname discrete.LR
#'@export
DLR <- function(raw.pvalues, pCDFlist, alpha = 0.05, zeta = 0.5, direction = "sd", critical.values = FALSE){
  return(discrete.LR(raw.pvalues, pCDFlist, alpha, zeta, direction, TRUE, critical.values))
}

#'@rdname discrete.LR
#'@export
NDLR <- function(raw.pvalues, pCDFlist, alpha = 0.05, zeta = 0.5, direction = "sd", critical.values = FALSE){
  return(discrete.LR(raw.pvalues, pCDFlist, alpha, zeta, direction, FALSE, critical.values))
}
