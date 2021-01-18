#'@title [HBR-\eqn{\lambda}] procedure
#'
#'@description
#'Apply the [HBR-\eqn{\lambda}] procedure,
#'with or without computing the critical constants,
#'to a set of p-values and their discrete support.
#'
#'@details
#'[DBR-lambda] is the discrete version of the [Blanchard-Roquain-lambda] procedure (see References),
#'the authors of the latter suggest to take \code{lambda = alpha} (see their Proposition 17),
#'which explains the choice of the default value here. 
#'
#'This version: 2019-06-18.
#'
#'@section References:
#'G. Blanchard and E. Roquain (2009). Adaptive false discovery rate control under independence and dependence. Journal of Machine Learning Research, 10, 2837-2871.
#'
#'@seealso
#'\code{\link{discrete.BH}}, \code{\link{DiscreteFDR}}
#'
#'@templateVar pvalues FALSE
#'@templateVar stepUp FALSE
#'@templateVar alpha TRUE
#'@templateVar sorted_pv FALSE
#'@templateVar support FALSE
#'@templateVar raw.pvalues TRUE
#'@templateVar pCDFlist TRUE
#'@templateVar direction FALSE
#'@templateVar ret.crit.consts TRUE
#'@templateVar lambda TRUE
#'@templateVar adaptive FALSE
#'@template param 
#'
#'@template example
#'@examples
#'
#'DBR.fast <- DBR(raw.pvalues, pCDFlist)
#'summary(DBR.fast)
#'DBR.crit <- DBR(raw.pvalues, pCDFlist, ret.crit.consts = TRUE)
#'summary(DBR.crit)
#'
#'@templateVar DBR TRUE
#'@template return
#'
#'@name DBR
#'@export
DBR <- function(raw.pvalues, pCDFlist, alpha = 0.05, lambda = NULL, ret.crit.consts = FALSE){
  # check arguments
  if(is.null(alpha) || is.na(alpha) || !is.numeric(alpha) || alpha < 0 || alpha > 1)
    stop("'alpha' must be a probability between 0 and 1!")
  if (is.null(lambda)){
    # if lambda is not provided, set lambda = alpha
    lambda <- alpha
  }else{
    if(is.na(lambda) || !is.numeric(lambda) || lambda < 0 || lambda > 1)
      stop("'lambda' must be a probability between 0 and 1!")
  }
  
  m <- length(raw.pvalues)
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
  #        Compute [DBR-lambda] significant p-values,
  #        their indices and the number of rejections
  #--------------------------------------------
  if(ret.crit.consts){
    # compute transformed support
    y <- kernel_DBR_crit(pCDFlist, pv.list.all, sorted.pvals, lambda, alpha)
    # find critical constants
    crit.constants <- y$crit.consts
    idx <- which(sorted.pvals <= crit.constants)   
  }
  else{
    # compute transformed p-values
    y <- kernel_DBR_fast(pCDFlist, sorted.pvals, lambda)
    idx <- which(y <= alpha)
  }
  m.rej <- length(idx)
  if(m.rej){
    idx <- which(pvec <= sorted.pvals[m.rej]) 
    pvec.rej <- raw.pvalues[idx]
  }else{
    idx <- integer(0)
    pvec.rej <- numeric(0)
  }
  #--------------------------------------------
  #       Create output list
  #--------------------------------------------
  output <- list(Rejected = pvec.rej, Indices = idx, Alpha = m.rej * alpha / m, Num.rejected = m.rej, Lambda = lambda)
  if(ret.crit.consts){
    # add critical values to output list
    output$Critical.values = crit.constants
    # compute adjusted p-values
    # recall that transformed p-values where max_i F_i(p) > lambda,
    # that is for indices > y$m.lambda, are set to 1
    pv.adj <- rev(cummin(rev(pmin(y$pval.transf, 1)))) # / c(seq_len(y$m.lambda), rep(1, m - y$m.lambda))
  }
  else{
    # compute adjusted p-values
    pv.adj <- rev(cummin(rev(pmin(y, 1))))
  }
  # add adjusted p-values to output list
  ro <- order(o)
  output$Adjusted = pv.adj[ro]
  
  # include details of the used algorithm as strings
  output$Method <- paste("Discrete Blanchard-Roquain procedure (lambda = ", lambda, ")", sep = "")
  output$Signif.level <- alpha
  output$Tuning <- lambda
  
  # original test data (often included, e.g. when using 'binom.test()')
  output$Data <- list()
  output$Data$raw.pvalues <- raw.pvalues
  output$Data$pCDFlist <- pCDFlist
  # object names of the data as strings
  output$Data$data.name <- paste(deparse(substitute(raw.pvalues)), "and", deparse(substitute(pCDFlist)))
  
  class(output) <- "DiscreteFDR"
  return(output)
  
  return(output)
}