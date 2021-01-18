#' SCORNET: A novel non-parametric survival curve estimator for the Electronic Health Record
#' 
#' Semi-Supervised Calibration of Risk with Noisy Event Times (SCORNET) is a consistent, non-parametric
#' survival curve estimator that boosts efficiency over existing non-parametric estimators
#' by (1) utilizing unlabeled patients in a semi-supervised fashion, and (2) leveraging
#' information-dense engineered EHR features to maximize unlabeled set imputation precision
#' See Ahuja et al. (2020) BioArxiv for details
#' @docType package
#' @name SCORNET-package
#' @keywords package
#' @useDynLib SCORNET
#' @import Rcpp
#' @import foreach
#' @import parallel
#' @import doParallel
#' @importFrom stats quasibinomial
#' @importFrom stats dnorm
#' @importFrom stats glm
#' @importFrom stats quantile
#' @importFrom stats sd
NULL
NULL
