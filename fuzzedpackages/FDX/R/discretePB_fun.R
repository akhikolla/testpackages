#'@name discrete.PB
#'
#'@title Discrete Poisson-Binomial procedure
#'
#'@description
#'Apply the [DPB] procedure, with or without computing the critical values, to
#'a set of p-values and their discrete support. A non-adaptive version is
#'available as well. Additionally, the user can choose between exact
#'computation of the Poisson-Binomial distribution or a refined normal
#'approximation.
#'
#'@details
#'\code{DPB} and \code{NDPB} are wrapper functions for \code{discrete.PB}.
#'The first one simply passes all its parameters to \code{discrete.PB} with
#'\code{adaptive = TRUE} and \code{NDPB} does the same with
#'\code{adaptive = FALSE}.
#'
#'@seealso
#'\code{\link{kernel}}, \code{\link{FDX-package}}, \code{\link{continuous.LR}},
#'\code{\link{continuous.GR}}, \code{\link{discrete.LR}}, 
#'\code{\link{discrete.GR}}, \code{\link{weighted.LR}},
#'\code{\link{weighted.GR}}, \code{\link{weighted.PB}}
#'
#'@templateVar raw.pvalues TRUE
#'@templateVar pCDFlist TRUE
#'@templateVar alpha TRUE
#'@templateVar zeta TRUE
#'@templateVar direction FALSE
#'@templateVar adaptive TRUE
#'@templateVar critical.values TRUE
#'@templateVar exact TRUE
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
#'DPB.fast <- DPB(raw.pvalues, pCDFlist)
#'summary(DPB.fast)
#'
#'DPB.crit <- DPB(raw.pvalues, pCDFlist, critical.values = TRUE)
#'summary(DPB.crit)
#'
#'NDPB.fast <- NDPB(raw.pvalues, pCDFlist)
#'summary(NDPB.fast)
#'
#'NDPB.crit <- NDPB(raw.pvalues, pCDFlist, critical.values = TRUE)
#'summary(NDPB.crit)
#'
#'@templateVar Critical.values TRUE
#'@template return
#'
#'@importFrom DiscreteFDR match.pvals
#'@export
discrete.PB <- function(raw.pvalues, pCDFlist, alpha = 0.05, zeta = 0.5, adaptive = TRUE, critical.values = FALSE, exact = TRUE){
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
  if(critical.values){
    # compute transformed support
    y <- kernel_DPB_crit(pCDFlist, pv.list.all, sorted.pvals, adaptive, alpha, zeta, exact)
    # find critical constants
    crit.constants <- y$crit.consts
    idx <- which(sorted.pvals > crit.constants)
  }else{
    # compute transformed sorted p-values
    y <- kernel_DPB_fast(pCDFlist, sorted.pvals, adaptive, alpha, exact)
    # determine significant (transformed) p-values
    idx <- which(y > zeta)
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
  #--------------------------------------------
  #       Create output S3 object
  #--------------------------------------------
  output <- list(Rejected = pvec.rej, Indices = idx, Num.rejected = m.rej)
  
  if(critical.values){
    y <- y$pval.transf
  }
  # compute adjusted p-values
  pv.adj = cummax(pmin(y, 1))
  # add adjusted p-values to output list
  ro <- order(o)
  output$Adjusted = pv.adj[ro]
  
  # add critical values to output list
  if(critical.values) output$Critical.values = crit.constants
  
  # include details of the used algorithm as strings
  alg <- "Discrete Poisson-Binomial procedure"
  alg <- if(exact) paste(alg, "(exact)") else paste(alg, "(normal approximation)")
  output$Method <- if(!adaptive) paste("Non-Adaptive", alg) else alg
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

#'@rdname discrete.PB
#'@export
DPB <- function(raw.pvalues, pCDFlist, alpha = 0.05, zeta = 0.5, critical.values = FALSE, exact = TRUE){
  return(discrete.PB(raw.pvalues, pCDFlist, alpha, zeta, TRUE, critical.values, exact))
}

#'@rdname discrete.PB
#'@export
NDPB <- function(raw.pvalues, pCDFlist, alpha = 0.05, zeta = 0.5, critical.values = FALSE, exact = TRUE){
  return(discrete.PB(raw.pvalues, pCDFlist, alpha, zeta, FALSE, critical.values, exact))
}