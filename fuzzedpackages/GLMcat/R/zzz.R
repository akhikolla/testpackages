## usethis namespace: start
#' @import Rcpp
# #' @export GLMadj
# #' @export GLMcum
# #' @export GLMseq
# #' @export GLMref
#' @export GLMcat
#' @export Discrete_CM
# #' @export predict_glmcat_Response
#' @export summary.glmcat
#' @export coef.glmcat
#' @export nobs_glmcat
#' @export logLik.glmcat
#' @export predict_glmcat
# #' @export ReferenceF
#' @useDynLib GLMcat, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @importFrom stats model.frame model.matrix as.formula pnorm printCoefmat
## usethis namespace: end

loadModule("GLMcatmodule", TRUE)
loadModule("discretemodule", TRUE)
# loadModule("discretemodule", TRUE)
# loadModule("cumulativemodule", TRUE)
# loadModule("exportmod", TRUE)
# loadModule("sequentialmodule", TRUE)
# loadModule("adjacentmodule", TRUE)
# loadModule("referencemodule", TRUE) # predict_glmcatION
