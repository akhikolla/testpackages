#' Generalized fiducial predictive distributions
#' @description Simulations of the generalized fiducial predictive 
#'   distributions.
#'
#' @param gfi a \code{\link{gfilmm}} object
#' @param newdata dataframe in which to look for variables with which to 
#'   predict, or \code{NULL} if the model is an intercept-only model without 
#'   random effect 
#'
#' @return A list with two fields: \code{FPD}, a dataframe containing the 
#'   simulations, and \code{WEIGHT}, their weight. This is a \code{gfilmm} 
#'   object.
#' 
#' @importFrom stats model.matrix rnorm
#' @importFrom utils head tail
#' @export
#' 
#' @note Actually the levels of the random effects given in \code{newdata} can 
#'   be different from the original levels. For instance, in the example 
#'   provided below, we enter \code{block = c("4","6")}, but we could also 
#'   enter \code{block = c("A","B")}, even though \code{"A"} and \code{"B"} 
#'   are not some levels of the \code{block} factor. Both options only mean 
#'   that the two observations to predict are in two different blocks.  
#'
#' @examples gfi <- gfilmm(
#'   ~ cbind(yield-0.1, yield+0.1), ~ N, ~ block, npk, 2000, nthreads = 2
#' )
#' fpd <- gfilmmPredictive(gfi, data.frame(N = c("0","1"), block = c("4","6")))
#' gfiSummary(fpd)
gfilmmPredictive <- function(gfi, newdata){
  if(is.null(newdata) || missing(newdata)){
    newdata <- as.data.frame(matrix(nrow = 1L, ncol = 0L))
  }
  if(anyDuplicated(newdata)){
    stop(
      "There are some duplicated rows in `newdata`."
    )
  }
  cvrts <- attr(gfi, "covariates")
  factors <- names(cvrts[["categorical"]])
  continuous <- names(cvrts[["continuous"]])
  nms   <- c(continuous, factors)
  if(any(!is.element(nms, names(newdata)))){
    stop(
      "`newdata` does not contain all the necessary variables."
    )
  }
  newdata <- droplevels(newdata)
  for(fact in factors){
    newdata[[fact]] <- as.factor(newdata[[fact]]) -> f
    levs <- cvrts[["categorical"]][[fact]]
    if(any(!is.element(levels(f), levs))){
      stop(
        "Found a factor in `newdata` with an unrecognized level."
      )
    }
    levels(newdata[[fact]]) <- levs
  }
  for(x in continuous){
    if(!is.numeric(newdata[[x]])){
      stop(
        sprintf("Variable `%s` should be numeric", x)
      )
    }
  }
  X <- model.matrix(attr(gfi, "fixed"), data = newdata) 
  RE2 <- getRE2(newdata, attr(gfi, "random"), check = FALSE)
  Z   <- getZ(RE2)
  N <- length(gfi[["WEIGHT"]])
  vertices <- t(as.matrix(gfi[["VERTEX"]]))
  neffects <- attr(gfi, "effects")
  fparams  <- head(vertices, neffects[["fixed"]])
  vars     <- tail(vertices, neffects[["random"]])
  E   <- vapply(RE2, nlevels, integer(1L))
  gauss <- matrix(rnorm(N*nrow(newdata)), nrow = N, ncol = nrow(newdata))
  out <- matrix(NA_real_, nrow = N, ncol = nrow(newdata))
  colnames(out) <- paste0("y", seq_len(nrow(newdata)))
  for(i in 1L:N){
    cholSigma <- chol(Z %*% (rep(vars[, i]^2, times = E) * t(Z)))
    Mu <- X %*% fparams[, i]
    out[i,] <- t(Mu) + t(gauss[i,]) %*% cholSigma # or gauss[i,,drop=FALSE] instead of t() => TODO: benchmarks
  }
  out <- list(FPD = as.data.frame(out), WEIGHT = gfi[["WEIGHT"]])
  class(out) <- c("gfilmm", "gfilmm.pred")
  out
}