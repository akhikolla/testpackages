#'@title [HSU], [HSD], [AHSU] and [AHSD] procedures
#'
#'@description
#'Apply the [HSU], [HSD], [AHSU] and [AHSD] procedures,
#'with or without computing the critical constants,
#'to a set of p-values and their discrete support.
#'
#'@details
#'\code{DBH} and \code{ADBH} are wrapper functions for \code{discrete.BH}. 
#'\code{DBH} simply passes all its parameters to \code{discrete.BH} with \code{adaptive = FALSE}. 
#'\code{ADBH} does the same with \code{adaptive = TRUE}.
#'
#'This version: 2019-06-18.
#'
#'@seealso
#'\code{\link{kernel}}, \code{\link{DiscreteFDR}}, \code{\link{DBR}}
#'
#'@templateVar pvalues FALSE
#'@templateVar stepUp FALSE
#'@templateVar alpha TRUE
#'@templateVar sorted_pv FALSE
#'@templateVar support FALSE
#'@templateVar raw.pvalues TRUE
#'@templateVar pCDFlist TRUE
#'@templateVar direction TRUE
#'@templateVar ret.crit.consts TRUE
#'@templateVar lambda FALSE
#'@templateVar adaptive TRUE
#'@template param 
#'
#'@template example
#'@examples
#'
#'DBH.su.fast <- DBH(raw.pvalues, pCDFlist)
#'summary(DBH.su.fast)
#'DBH.sd.fast <- DBH(raw.pvalues, pCDFlist, direction = "sd")
#'DBH.sd.fast$Adjusted
#'summary(DBH.sd.fast)
#'
#'DBH.su.crit <- DBH(raw.pvalues, pCDFlist, ret.crit.consts = TRUE)
#'summary(DBH.su.crit)
#'DBH.sd.crit <- DBH(raw.pvalues, pCDFlist, direction = "sd", ret.crit.consts = TRUE)
#'DBH.sd.crit$Adjusted
#'summary(DBH.sd.crit)
#'
#'ADBH.su.fast <- ADBH(raw.pvalues, pCDFlist)
#'summary(ADBH.su.fast)
#'ADBH.sd.fast <- ADBH(raw.pvalues, pCDFlist, direction = "sd")
#'ADBH.sd.fast$Adjusted
#'summary(ADBH.sd.fast)
#'
#'ADBH.su.crit <- ADBH(raw.pvalues, pCDFlist, ret.crit.consts = TRUE)
#'summary(ADBH.su.crit)
#'ADBH.sd.crit <- ADBH(raw.pvalues, pCDFlist, direction = "sd", ret.crit.consts = TRUE)
#'ADBH.sd.crit$Adjusted
#'summary(ADBH.sd.crit)
#'
#'@templateVar DBR FALSE
#'@template return
#'
#'@name discrete.BH
NULL

#'@rdname discrete.BH
#'@export
discrete.BH <- function(raw.pvalues, pCDFlist, alpha = 0.05, direction = "su", adaptive = FALSE, ret.crit.consts = FALSE){
  # check arguments
  if(is.null(alpha) || is.na(alpha) || !is.numeric(alpha) || alpha < 0 || alpha > 1)
    stop("'alpha' must be a probability between 0 and 1!")
  
  m <- length(raw.pvalues)
  if(m != length(pCDFlist)) stop("The lengths of 'raw.pvalues' and 'pCDFlist' must be equal!")
  
  direction <- match.arg(direction, c("su", "sd"))
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
    if(ret.crit.consts){
      if(adaptive){
        # compute transformed support
        y <- kernel_ADBH_crit(pCDFlist, pv.list.all, sorted.pvals, TRUE, alpha)
      }
      else{
        # compute transformed support
        y <- kernel_DBH_crit(pCDFlist, pv.list.all, sorted.pvals, TRUE, alpha)
      }
      # find critical constants
      crit.constants <- y$crit.consts
      idx <- which(sorted.pvals <= crit.constants)
    }
    else{
      if(adaptive){
        # compute transformed observed p-values
        y <- kernel_ADBH_fast(pCDFlist, sorted.pvals, TRUE, alpha, pv.list.all)
      }
      else{
        # compute transformed observed p-values
        y <- kernel_DBH_fast(pCDFlist, sorted.pvals, TRUE, alpha, pv.list.all)
      }
      # determine significant (transformed) p-values
      if(length(y)){
        idx <- which(y <= 1:length(y) * alpha)
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
    if(ret.crit.consts){
      if(adaptive){
        # compute transformed support
        y <- kernel_ADBH_crit(pCDFlist, pv.list.all, sorted.pvals, FALSE, alpha)
      }
      else{
        # compute transformed support
        y <- kernel_DBH_crit(pCDFlist, pv.list.all, sorted.pvals, FALSE, alpha)
      }
      # find critical constants
      crit.constants <- y$crit.consts
      idx <- which(sorted.pvals > crit.constants)
    }
    else{
      if(adaptive){
        # compute transformed sorted p-values
        y <- kernel_ADBH_fast(pCDFlist, sorted.pvals, FALSE)
      }
      else{
        # compute transformed sorted p-values
        y <- kernel_DBH_fast(pCDFlist, sorted.pvals, FALSE)
      }
      # determine significant (transformed) p-values
      idx <- which(y > 1:m * alpha) 
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
    if(ret.crit.consts){
      y <- y$pval.transf
    }
    # compute adjusted p-values
    pv.adj = cummax(pmin(y / 1:m, 1))
    # add adjusted p-values to output list
    ro <- order(o)
    output$Adjusted = pv.adj[ro]
  }
  # add critical values to output list
  if(ret.crit.consts) output$Critical.values = crit.constants
  
  # include details of the used algorithm as strings
  alg <- "Discrete Benjamini-Hochberg procedure"
  alg <- if(adaptive) paste("Adaptive", alg) else alg
  output$Method <- paste(alg, switch(direction, su = "(step-up)", sd = "(step-down)"))
  output$Signif.level <- alpha
  
  # original test data (often included, e.g. when using 'binom.test()')
  output$Data <- list()
  output$Data$raw.pvalues <- raw.pvalues
  output$Data$pCDFlist <- pCDFlist
  # object names of the data as strings
  output$Data$data.name <- paste(deparse(substitute(raw.pvalues)), "and", deparse(substitute(pCDFlist)))
  
  class(output) <- "DiscreteFDR"
  return(output)
}

#'@rdname discrete.BH
#'@export
DBH <- function(raw.pvalues, pCDFlist, alpha = 0.05, direction = "su", ret.crit.consts = FALSE){
  return(discrete.BH(raw.pvalues, pCDFlist, alpha, direction, adaptive = FALSE, ret.crit.consts))
}

#'@rdname discrete.BH
#'@export
ADBH <- function(raw.pvalues, pCDFlist, alpha = 0.05, direction = "su", ret.crit.consts = FALSE){
  return(discrete.BH(raw.pvalues, pCDFlist, alpha, direction, adaptive = TRUE, ret.crit.consts))
}
