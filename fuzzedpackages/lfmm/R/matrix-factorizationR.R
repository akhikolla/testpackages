#' R package : Fast and Accurate statistical methods for adjusting confounding factors
#' in association studies.
#'
#'
#' @docType package
#'
#' @name lfmm
#' @description Implements statistical methods for adjusting confounding factors
#' in association studies.
#' @references 
#' Caye, K., B. Jumentier, J. Lepeule, and O. François, 2019  LFMM 2: fast and accurate inference of gene-environment associations in genome-widestudies. Mol. Biol. Evol. 36: 852–860.https://doi.org/10.1093/molbev/msz008
#' 
#' B. Jumentier, Caye, K., J. Lepeule, and O. François, 2019 Sparse latent factor regression models for genome-wide and epigenome-wide association studies (in prep)
#' @importFrom Rcpp evalCpp
#' @importFrom foreach foreach %:% %do% %dopar%
#' @useDynLib lfmm
#' @import RcppEigen
NULL

# Fit the model
# @export
lfmm_fit <- function(m, dat, ...) {
  UseMethod("lfmm_fit")
}

# Fit the model when latent factor loadings are known
# @export
lfmm_fit_knowing_loadings <- function(m, dat, ...) {
  UseMethod("lfmm_fit_knowing_loadings")
}

# Cross validation
# @export
# lfmm_CV <- function(m, dat, n.fold.row, n.fold.col, ...) {
#   UseMethod("lfmm_CV")
# }

# Impute Y with a fitted model.
# @export
lfmm_impute <- function(m, dat, ...) {
  UseMethod("lfmm_impute")
}

# Compute the residual error
# @export
lfmm_residual_error2 <- function(m, dat, ...) {
  UseMethod("lfmm_residual_error2")
}
