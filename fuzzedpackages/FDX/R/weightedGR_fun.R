#'@name weighted.GR
#'
#'@title Weighted Guo-Romano Procedure
#'
#'@description
#'Apply the weighted [wGR] procedure, with or without computing the
#'critical values, to a set of p-values. Both arithmetic and geometric
#'weighting are available.
#'
#'@details
#'\code{wGR.AM} and \code{wGR.GM} are wrapper functions for \code{weighted.GR}.
#'The first one simply passes all its parameters to \code{weighted.GR} with
#'\code{weighting.method = "AM"} and \code{wGR.GM} does the same with
#'\code{weighting.method = "GM"}.
#'
#'@seealso
#'\code{\link{kernel}}, \code{\link{FDX-package}}, \code{\link{continuous.LR}},
#'\code{\link{continuous.GR}}, \code{\link{discrete.LR}}, 
#'\code{\link{discrete.GR}}, \code{\link{discrete.PB}}, 
#'\code{\link{weighted.LR}}, \code{\link{weighted.PB}}
#'
#'@templateVar raw.pvalues TRUE
#'@templateVar pCDFlist FALSE
#'@templateVar alpha TRUE
#'@templateVar zeta TRUE
#'@templateVar direction FALSE
#'@templateVar adaptive FALSE
#'@templateVar critical.values TRUE
#'@templateVar exact FALSE
#'@templateVar pvalues FALSE
#'@templateVar sorted_pv FALSE
#'@templateVar stepUp FALSE
#'@templateVar support FALSE
#'@templateVar weights TRUE
#'@templateVar weighting.method TRUE
#'@template param
#' 
#'@section References:
#' S. Döhler and E. Roquain (2019). Controlling False Discovery Exceedance for
#' Heterogeneous Tests.
#' \href{https://arxiv.org/abs/1912.04607v1}{arXiv:1912.04607v1}.
#'
#'@template exampleWeighted
#'@examples
#'
#'wGR.AM.fast <- wGR.AM(raw.pvalues.weighted, weights)
#'summary(wGR.AM.fast)
#'
#'wGR.AM.crit <- wGR.AM(raw.pvalues.weighted, weights, critical.values = TRUE)
#'summary(wGR.AM.crit)
#'
#'wGR.GM.fast <- wGR.GM(raw.pvalues.weighted, weights)
#'summary(wGR.GM.fast)
#'
#'wGR.GM.crit <- wGR.GM(raw.pvalues.weighted, weights, critical.values = TRUE)
#'summary(wGR.GM.crit)
#'
#'@templateVar Critical.values TRUE
#'@templateVar Adaptive FALSE
#'@templateVar Weighting TRUE
#'@template return
#'
#'@importFrom pracma fzero
#'@importFrom stats pbinom qbeta
#'@export
weighted.GR <- function(raw.pvalues, weights, alpha = 0.05, zeta = 0.5, weighting.method = "AM", critical.values = FALSE){
  #--------------------------------------------
  #       check arguments
  #--------------------------------------------
  if(is.null(raw.pvalues) || !is.numeric(raw.pvalues) || all(is.na(raw.pvalues)) || any(raw.pvalues < 0 | raw.pvalues > 1))
    stop("All values of 'raw.pvalues' must be probabilities between 0 and 1!")
  
  m <- length(raw.pvalues)
  if(m != length(weights)) stop("The lengths of 'raw.pvalues' and 'weights' must be equal!")
  
  if(is.null(weights) || !is.numeric(weights) || all(is.na(weights)))
    stop("All values of 'weights' must be numeric!")
  
  if(is.null(alpha) || is.na(alpha) || !is.numeric(alpha) || alpha < 0 || alpha > 1)
    stop("'alpha' must be a probability between 0 and 1!")
  
  if(is.null(zeta) || is.na(zeta) || !is.numeric(zeta) || zeta < 0 || zeta > 1)
    stop("'zeta' must be a probability between 0 and 1!")
  
  weighting.method <- match.arg(weighting.method, c("AM", "GM"))
  #--------------------------------------------
  #       remove NA's from raw p.values
  #--------------------------------------------
  non.missing <- !is.na(raw.pvalues) & !is.na(weights)
  raw.pvalues.ok <- raw.pvalues[non.missing]
  weights.ok <- weights[non.missing]
  #--------------------------------------------
  #       number of tests, m(l) and [alpha * m] + 1
  #--------------------------------------------
  m <- length(raw.pvalues.ok)
  a <- floor(alpha * 1:m)
  m.l <- m - 1:m + a + 1
  #--------------------------------------------
  #       Rescale weights
  #--------------------------------------------
  weights.rescaled <- weights.ok / mean(weights.ok)
  weights.decreasing <- sort(weights.rescaled, decreasing = TRUE)
  #--------------------------------------------
  #       Compute weighted p-values
  #--------------------------------------------
  switch(weighting.method,
         AM = {qvalues <- pmin(raw.pvalues.ok / weights.rescaled, 1)},
         GM = {qvalues <- geom_weight(raw.pvalues.ok, 1 / weights.rescaled)}
  )
  #--------------------------------------------
  #       Determine sort order and do sorting
  #--------------------------------------------
  o <- order(qvalues)
  ro <- order(o)
  qvalues.sorted <- qvalues[o]
  #--------------------------------------------
  #        Compute [wLR-AM] or [wLR-GM] significant p-values,
  #        their indices and the number of rejections
  #--------------------------------------------
  # Compute transformed weighted p-values
  y <- pmin(1, cummax(kernel_wGR_fast(qvalues.sorted, weights.decreasing, alpha, weighting.method == "GM")))
  # determine significant (transformed) p-values
  idx <- which(y > zeta)
  if(length(idx)){
    m.rej <- min(idx) - 1
    if (m.rej){
      # determine significant (observed) p-values in sorted.pvals
      idx <- which(qvalues <= qvalues.sorted[m.rej])
      pvec.rej <- raw.pvalues.ok[idx]
    }else{
      idx <- numeric(0)
      pvec.rej <- numeric(0)
    }
  }else{
    m.rej <- m
    idx <- 1:m
    pvec.rej <- raw.pvalues.ok
  }
  # find critical values
  if(critical.values){
    if(weighting.method == "AM"){
      crit <- numeric(m)
      for(l in 1:m) crit[l] = fzero(function(x, w, b, z){ppbinom(b, pmin(w * x, 1), method = "GeoMeanCounter", lower.tail = FALSE) - z}, c(max(crit), 1), w = weights.decreasing[1:m.l[l]], b = a[l], z = zeta, tol = .Machine$double.neg.eps)$x
    }else{
      crit <- geom_weight(qbeta(zeta, a + 1, m.l - a), (1:m / cumsum(weights.decreasing))[m.l])
    }
  }
  
  #--------------------------------------------
  #       Create output S3 object
  #--------------------------------------------
  output <- list(Rejected = pvec.rej, Indices = idx, Num.rejected = m.rej, Adjusted = y[ro], Weighted = qvalues)
  
  # add critical values to output list
  if(critical.values) output$Critical.values <- crit
  
  # include details of the used algorithm as strings
  alg <- "Weighted Guo-Romano procedure"
  output$Method <- if(weighting.method == "AM") paste(alg, "(arithmetic weighting)") else paste(alg, "(geometric weighting)")
  output$FDP.threshold <- alpha
  output$Exceedance.probability <- zeta
  output$Weighting <- ifelse(weighting.method == "AM", "Arithmetic", "Geometric")
  
  # original test data
  output$Data <- list()
  output$Data$raw.pvalues <- raw.pvalues
  output$Data$weights <- weights
  # object names of the data as strings
  output$Data$data.name <- deparse(substitute(raw.pvalues))
  
  class(output) <- "FDX"
  return(output)
}

#'@rdname weighted.GR
#'@export
wGR.AM <- function(raw.pvalues, weights, alpha = 0.05, zeta = 0.5, critical.values = FALSE){
  return(weighted.GR(raw.pvalues, weights, alpha, zeta, "AM", critical.values))
}

#'@rdname weighted.GR
#'@export
wGR.GM <- function(raw.pvalues, weights, alpha = 0.05, zeta = 0.5, critical.values = FALSE){
  return(weighted.GR(raw.pvalues, weights, alpha, zeta, "GM", critical.values))
}