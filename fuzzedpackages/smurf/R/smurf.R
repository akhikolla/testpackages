###############################################
#
# smurf package
#
###############################################


#############################
# Import functions from packages

#' @import catdata
#'
#' @importFrom glmnet glmnet
#' 
#' @importFrom graphics abline
#' @importFrom graphics plot
#' @importFrom graphics points
#' @importFrom graphics segments
#' 
#' @importFrom MASS ginv
#' 
#' @importFrom Matrix colSums
#' @importFrom Matrix isSymmetric
#' @importFrom Matrix Matrix
#' @importFrom Matrix rowSums
#' @importFrom Matrix rankMatrix
#' @importFrom Matrix sparse.model.matrix
#' @importFrom Matrix triu
#' 
#' @importFrom methods as
#' 
#' @importFrom mgcv gam
#' @importFrom mgcv predict.gam
#' 
#' @importFrom parallel clusterExport
#' @importFrom parallel detectCores
#' @importFrom parallel makeCluster
#' @importFrom parallel parLapply
#' @importFrom parallel stopCluster
#' 
# Import function from Rcpp to import registered functions
#' @importFrom Rcpp evalCpp
#' 
#' @importFrom RColorBrewer brewer.pal
#' 
#' @importFrom speedglm speedglm.wfit
#' 
#' @importFrom stats as.formula
#' @importFrom stats contrasts
#' @importFrom stats coef
#' @importFrom stats coefficients 
#' @importFrom stats delete.response
#' @importFrom stats fitted
#' @importFrom stats glm
#' @importFrom stats glm.fit
#' @importFrom stats median
#' @importFrom stats model.frame
#' @importFrom stats model.offset
#' @importFrom stats model.response
#' @importFrom stats optim
#' @importFrom stats predict
#' @importFrom stats relevel 
#' @importFrom stats resid
#' @importFrom stats residuals
#' @importFrom stats sd
#' @importFrom stats terms
#' @importFrom stats weighted.mean

#############################
# Numerical tolerance
eps_num <- min(sqrt(.Machine$double.eps), 100 * .Machine$double.eps)

#############################
# Make C code recognisable
#' @useDynLib smurf, .registration = TRUE 
"_PACKAGE"

#############################
# Unload package DLL when package is unloaded
.onUnload <- function (libpath) {
  library.dynam.unload("smurf", libpath)
}
