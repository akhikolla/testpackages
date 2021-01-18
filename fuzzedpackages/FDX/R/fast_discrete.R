#'@title Fast application of discrete procedures
#'@name fast.Discrete
#'
#'@importFrom DiscreteFDR fisher.pvalues.support
#'
#'@description
#'Applies the [DLR], [DGR] or [DPB] procedures, without computing the critical
#'values, to a data set of 2 x 2 contingency tables using Fisher's exact test.
#'
#'
#'@param counts        a data frame of 2 or 4 columns and any number of lines,
#'                     each line representing a 2 x 2 contingency table to
#'                     test. The number of columns and what they must contain
#'                     depend on the value of the \code{input} argument, see
#'                     Details of \code{\link{fisher.pvalues.support}}.
#'@param alternative   same argument as in \code{\link{fisher.test}}. The three
#'                     possible values are \code{"greater"} (default),
#'                     \code{"two.sided"} or \code{"less"}; may be abbreviated.
#'@param input         the format of the input data frame, see Details of
#'                     \code{\link[DiscreteFDR]{fisher.pvalues.support}}. The
#'                     three possible values are \code{"noassoc"} (default),
#'                     \code{"marginal"} or \code{"HG2011"}; may be 
#'                     abbreviated.
#'
#'@templateVar raw.pvalues FALSE
#'@templateVar pCDFlist FALSE
#'@templateVar alpha TRUE
#'@templateVar zeta TRUE
#'@templateVar direction TRUE
#'@templateVar adaptive TRUE
#'@templateVar critical.values FALSE
#'@templateVar exact TRUE
#'@templateVar pvalues FALSE
#'@templateVar sorted_pv FALSE
#'@templateVar stepUp FALSE
#'@templateVar support FALSE
#'@templateVar weights FALSE
#'@templateVar weighting.method FALSE
#'@template param
#'
#'@examples
#'
#'X1 <- c(4, 2, 2, 14, 6, 9, 4, 0, 1)
#'X2 <- c(0, 0, 1, 3, 2, 1, 2, 2, 2)
#'N1 <- rep(148, 9)
#'N2 <- rep(132, 9)
#'Y1 <- N1 - X1
#'Y2 <- N2 - X2
#'df <- data.frame(X1, Y1, X2, Y2)
#'df
#'
#'DLR.sd <- fast.Discrete.LR(counts = df, input = "noassoc")
#'DLR.sd$Adjusted
#'summary(DLR.sd)
#'DLR.su <- fast.Discrete.LR(counts = df, input = "noassoc", direction = "su")
#'summary(DLR.su)
#'
#'NDLR.sd <- fast.Discrete.LR(counts = df, input = "noassoc", adaptive = FALSE)
#'NDLR.sd$Adjusted
#'summary(NDLR.sd)
#'NDLR.su <- fast.Discrete.LR(counts = df, input = "noassoc", direction = "su", adaptive = FALSE)
#'summary(NDLR.su)
#'
#'DGR <- fast.Discrete.GR(counts = df, input = "noassoc")
#'DGR$Adjusted
#'summary(DGR)
#'
#'NDGR <- fast.Discrete.GR(counts = df, input = "noassoc", adaptive = FALSE)
#'NDGR$Adjusted
#'summary(NDGR)
#'
#'DPB <- fast.Discrete.PB(counts = df, input = "noassoc")
#'DPB$Adjusted
#'summary(DPB)
#'
#'NDPB <- fast.Discrete.PB(counts = df, input = "noassoc", adaptive = FALSE)
#'NDPB$Adjusted
#'summary(NDPB)
#'
#'@templateVar Critical.values FALSE
#'@templateVar Adaptive TRUE
#'@templateVar Weighting FALSE
#'@template return
#'
#'@export
fast.Discrete.LR <- function(counts, alternative = "greater", input = "noassoc", alpha = 0.05, zeta = 0.5, direction = "sd", adaptive = TRUE){
  data.formatted <- fisher.pvalues.support(counts, alternative, input)
  raw.pvalues <- data.formatted$raw
  pCDFlist <- data.formatted$support
  
  out <- discrete.LR(raw.pvalues, pCDFlist, alpha, zeta, direction, adaptive, FALSE)
  out$Data$data.name <- deparse(substitute(counts)) 
  
  return(out)
}

#'@rdname fast.Discrete
#'@export
fast.Discrete.PB <- function(counts, alternative = "greater", input = "noassoc", alpha = 0.05, zeta = 0.5, adaptive = TRUE, exact = FALSE){
  data.formatted <- fisher.pvalues.support(counts, alternative, input)
  raw.pvalues <- data.formatted$raw
  pCDFlist <- data.formatted$support
  
  out <- discrete.PB(raw.pvalues, pCDFlist, alpha, zeta, adaptive, FALSE, exact)
  out$Data$data.name <- deparse(substitute(counts)) 
  
  return(out)
}

#'@rdname fast.Discrete
#'@export
fast.Discrete.GR <- function(counts, alternative = "greater", input = "noassoc", alpha = 0.05, zeta = 0.5, adaptive = TRUE){
  data.formatted <- fisher.pvalues.support(counts, alternative, input)
  raw.pvalues <- data.formatted$raw
  pCDFlist <- data.formatted$support
  
  out <- discrete.GR(raw.pvalues, pCDFlist, alpha, zeta, adaptive, FALSE)
  out$Data$data.name <- deparse(substitute(counts)) 
  
  return(out)
}