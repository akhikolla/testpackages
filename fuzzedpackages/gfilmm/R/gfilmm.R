#' Generalized fiducial inference
#' @description Samples the fiducial distributions.
#'
#' @param y a right-sided formula of the form \code{~ cbind(lower,upper)} for 
#'   the interval data 
#' @param fixed a right-sided formula for the fixed effects
#' @param random a right-sided formula for the random effects, or \code{NULL} 
#'   for no random effect
#' @param data the data, a dataframe
#' @param N desired number of simulations
#' @param thresh threshold, default \code{N/2}; for experts only
#' @param long logical, whether to use long doubles instead of doubles in the 
#'   algorithm
#' @param seed the seed for the C++ random numbers generator, a positive 
#'   integer, or \code{NULL} to use a random seed 
#' @param nthreads number of threads to run the algorithm with parallelized 
#'   blocks of code
#' @param x a \code{gfilmm} object
#' @param ... ignored
#'
#' @return A list with two components: a dataframe \code{VERTEX}, and a vector 
#'   \code{WEIGHT}. It has class \code{gfilmm}.
#' 
#' @importFrom stats model.matrix terms.formula as.formula
#' @importFrom utils head
#' @importFrom rgr gx.sort.df
#' @importFrom Matrix sparseMatrix
#' @importFrom parallel detectCores
#' @export
#' 
#' @references J. Cisewski and J.Hannig. 
#'   \emph{Generalized fiducial inference for normal linear mixed models}. 
#'   The Annals of Statistics 2012, Vol. 40, No. 4, 2102–2127.
#'
#' @examples h <- 0.01
#' gfi <- gfilmm(
#'   ~ cbind(yield-h, yield+h), ~ 1, ~ block, data = npk, N = 5000, nthreads = 2
#' )
#' # fiducial cumulative distribution function of the intercept:
#' Fintercept <- gfiCDF(~ `(Intercept)`, gfi)
#' plot(Fintercept, xlim = c(40, 65))
#' # fiducial confidence interval of the intercept:
#' gfiConfInt(~ `(Intercept)`, gfi)
#' # fiducial density function of the intercept:
#' library(kde1d)
#' kfit <- kde1d(gfi$VERTEX[["(Intercept)"]], weights = gfi$WEIGHT)
#' curve(dkde1d(x, kfit), from = 45, to = 65)
gfilmm <- function(
  y, fixed, random, data, 
  N, thresh = N/2, long = FALSE, seed = NULL, 
  nthreads = parallel::detectCores()
){
  if(inSolaris()) nthreads <- 1L
  seed <- if(is.null(seed)){
    sample.int(1000000L, 1L)
  }else{
    as.integer(abs(seed))
  }
  data <- droplevels(data)
  if(!is.null(random)){
    factors <- rownames(attr(terms.formula(random), "factors"))
    frmla <- as.formula(paste0("~ ", paste0(factors, collapse = " + ")))
    data <- gx.sort.df(frmla, data)
  }
  Y <- f_eval_rhs(y, data = data)
  if(!is.matrix(Y) || ncol(Y) != 2L){
    stop(
      "`y` must be a right-sided formula of the form `~ cbind(lower, upper)`."
    )
  }
  if(!is.numeric(Y)){
    stop(
      "Invalid `y` argument: the response must be given by two numerical ",
      "columns of the data."
    )
  }
  n <- nrow(Y)
  yl <- Y[,1L]; yu <- Y[,2L]
  if(any(yl >= yu)){
    stop(
      "Invalid interval data: found some values in the first column higher ",
      "than the corresponding values in the second column."
    )
  }
  FE <- model.matrix(fixed, data = data)
  fullrank <- qr(FE)$rank == ncol(FE)
  if(!fullrank){
    stop(
      "The design matrix of the fixed effects is not of full rank."
    )
  }
  RE2 <- getRE2(data, random, check = TRUE)
  tlabs <- head(names(RE2), -1L)
  Z <- getZ(RE2) 
  ij <- which(Z == 1L, arr.ind = TRUE)
  Z <- sparseMatrix(i = ij[,1L], j = ij[,2L], dims = dim(Z), x = 1.0)
  E <- vapply(RE2, nlevels, integer(1L))
  RE2 <- vapply(RE2, recode, integer(n))#vapply(RE2, as.integer, integer(n)) - 1L # recode
  gfi <- if(long){
    gfilmm_long(yl, yu, FE, Z, RE2, E, N, thresh, seed, nthreads)
  }else{
    gfilmm_double(yl, yu, FE, Z, RE2, E, N, thresh, seed, nthreads)
  }
  rownames(gfi[["VERTEX"]]) <-
    c(colnames(FE), paste0("sigma_", c(tlabs, "error")))
  gfi[["VERTEX"]] <- as.data.frame(t(gfi[["VERTEX"]]))
  attr(gfi, "effects") <- c(fixed = ncol(FE), random = ncol(RE2))
  attr(gfi, "covariates") <- getCovariates(data, fixed)
  attr(gfi, "fixed") <- fixed
  attr(gfi, "random") <- random
  class(gfi) <- "gfilmm"
  gfi
}

#' @rdname gfilmm
#' @importFrom stats setNames
#' @importFrom utils capture.output str
#' @export
print.gfilmm <- function(x, ...){
  nms <- names(x)
  attributes(x) <- NULL
  cat(
    "`gfilmm` object", "\n",
    "---------------", "\n",
    sep = ""
  )
  cat(capture.output(str(setNames(x, nms))), sep = "\n")
}

