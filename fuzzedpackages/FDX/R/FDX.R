#'@name FDX-package
#'
#'@docType package
#'@import Rcpp
#'@importFrom Rcpp evalCpp sourceCpp
#'@useDynLib FDX
#'@title False Discovery Exceedance (FDX) Control for Heterogeneous and Discrete Tests
#'
#'@description
#'This package implements the [HLR], [HGR] and [HPB] procedures for both
#'heterogeneous and discrete tests (see Reference). 
#'
#'@details
#'The functions are reorganized from the reference paper in the following way.
#'\code{discrete.LR} (for Discrete Lehmann-Romano) implements [DLR],
#'\code{discrete.GR} (for Discrete Guo-Romano) implements [DGR] and
#'\code{discrete.PB} (for Discrete Poisson-Binomial) implements [DPB].
#'\code{DLR} and \code{NDLR} are wrappers for \code{discrete.LR} to access
#'[DLR] and its non-adaptive version directly. Likewise, \code{DGR},
#'\code{NDGR}, \code{DPB} and \code{NDPB} are wrappers for
#'\code{discrete.GR} and \code{discrete.PB}, respectively. Their main
#'parameters are a vector of raw observed p-values and a list of the same
#'length, whose elements are the discrete supports of the CDFs of the p-values.
#'
#'In the same fashion, \code{weighted.LR} (for Weighted Lehmann-Romano),
#'\code{weighted.GR} (for Weighted Guo-Romano) and \code{weighted.PB}
#'(for Weighted Poisson-Binomial) implement [wLR], [wGR] and [wGR],
#'respectively. They also possess wrapper functions, namely \code{wLR.AM},
#'\code{wGR.AM} and \code{wPB.AM} for arithmetic weighting, and \code{wLR.GM},
#'\code{wPB.GM} and \code{wPB.GM} for geometric weighting.
#' 
#'The functions \code{fast.Discrete.LR}, \code{fast.Discrete.GR}
#'and \code{fast.Discrete.PB} are wrappers for
#'\code{\link[DiscreteFDR]{fisher.pvalues.support}} and \code{discrete.LR},
#'\code{discrete.GR} and \code{discrete.PB}, respectively, which allow to apply
#'discrete procedures directly to a data set of contingency tables.
#' 
#'@section References:
#' S. Döhler and E. Roquain (2019). Controlling False Discovery Exceedance for
#' Heterogeneous Tests.
#' \href{https://arxiv.org/abs/1912.04607v1}{arXiv:1912.04607v1}.
NULL
