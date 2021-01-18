#' @title Calculates bivariate Pearson, Spearman, polychoric, and polyserial correlation coefficients
#'
#' @description Calculates bivariate Pearson, Spearman, polychoric, and polyserial correlation
#' coefficients in weighted or unweighted form, on discrete or continuous variables. Also 
#' calculates tetrachoric and biserial correlation coefficients as described below.
#'
#' @param x          a numeric (or numeric factor in case of polychoric) vector or an object that can be
#'                   coerced to a numeric or factor vector.
#' @param y          a numeric vector (or factor in case of polychoric and polyserial) or an object that
#'                   can be coerced to a numeric or factor vector.
#' @param method     a character string indicating which correlation coefficient is
#'                   to be computed. These include "Pearson" (default), "Spearman", "Polychoric", or "Polyserial".
#'                   For tetrachoric use "Polychoric" and for biserial use "Polyserial".
#' @param weights    a numeric vector of weights. By default, the unweighted correlation coefficient is calculated
#'                   by setting the weights to a vector of all 1s.
#' @param ML         a Boolean value indicating if full Maximum Likelihood (ML) is to be used (polyserial and polychoric only,
#'                   has no effect on Pearson or Spearman results). This substantially increases the
#'                   compute time. See the 'wCorr Arguments' vignette for a description of the effect of this argument.
#' @param fast       a Boolean value indicating if the Rcpp methods should be used. Setting this value to FALSE
#'                   uses the pure R implementation and is included primarily for comparing the implementations
#'                   to each other. See the 'wCorr Arguments' vignette for a description of the effect of this argument.
#'
#' @details 
#' In case of polyserial, x must be the observed ordinal variable, and y the observed continuous variable. For
#' polychoric, both must be categorical. The correlation methods are calculated as described in the 'wCorr Formulas'
#' vignette.
#'
#' For Spearman the data is first ranked and then a Pearson type correlation coefficient is calculated on
#' the result. The ranking method gives averages for ties.
#'
#' The details of computation are given in the 'wCorr Formulas' vignette.
#'
#' @return
#' A scalar that is the estimated correlation.
#'
#' @references
#'  Polyserial computation based on the likelihood function in Cox, N. R. (1974), "Estimation of the Correlation between a Continuous and a Discrete Variable." Biometrics, 30 (1), pp 171-178.
#'
#' Polychoric computation based on the likelihood function in Olsson, U. (1979) "Maximum Likelihood Estimation of the Polychoric Correlation Coefficient." Psyhometrika, 44 (4), pp 443-460.
#' 
#' The weighted Pearson formula appears in many places, including the "correlate" function in Stata Corp, Stata Statistical Software: Release 8. College Station, TX: Stata Corp LP, 2003.
#' 
#'
#' @examples
#' # run a polyserial correlation
#' attach(mtcars)
#' weightedCorr(gear, x=cyl, method="polyserial")
#' # weight by MPG
#' weightedCorr(y=gear, x=cyl, method="polyserial", weights=mpg)
#' # unweight
#' weightedCorr(y=gear, x=cyl, method="polyserial")
#'
#' # run a polychoric correlation
#' weightedCorr(gear, x=cyl, method="polychoric")
#' # weight by MPG
#' weightedCorr(y=gear, x=cyl, method="polychoric", weights=mpg)
#' # unwiehgted
#' weightedCorr(y=gear, x=cyl, method="polychoric")
#' detach(mtcars)
#'
#' @seealso \ifelse{latex}{\code{cor}}{\code{\link[stats]{cor}}}
#'
#' @export
weightedCorr <- function(x, y, method = c("Pearson", "Spearman", "Polyserial", "Polychoric"), weights=rep(1,length(x)), ML=FALSE, fast=TRUE) {
  x <- as.numeric(x)
  y <- as.numeric(y)
  weights <- as.numeric(weights)
  if(!is.vector(x)) stop(paste0("The argument ",sQuote("x"), " must be a vector."))
  if(!is.vector(y)) stop(paste0("The argument ",sQuote("y"), " must be a vector."))
  if(!is.vector(weights)) stop(paste0("The argument ",sQuote("weights"), " must be a vector."))
  if(length(x) != length(y)) stop(paste0("The vectors ", sQuote("x"), ", ", sQuote("y"), ", must be the same length."))
  if(length(x) != length(weights)) stop(paste0("The vectors ", sQuote("x"), ", ", sQuote("y"), ", and ", sQuote("weights") ," must all be of the same length."))

  value <- 0
  foundMethod <- FALSE
  #if (method == "Polyserial") {
  #  value <- polys(x, y, weights)
  #}
  method0 <- method
  method <- tolower(method)

  if(method == "polyserial") {
    if(length(unique(y)) == length(y) & length(unique(x)) < length(x)) {
      stop(paste0("Check argument definitions for ", sQuote("y"), " and ", sQuote("x") ,". The number of levels in the discrete variable ",sQuote("y")," is equal to the number of observations while the number of levels in continuous variable ",  sQuote("x")," is less than the number of observations. Try transposing these two arguments."))
    }
    if(length(unique(y)) > length(unique(x))) {
      warning(paste0("Check argument definitions for ", sQuote("y"), " and ", sQuote("x") ,". The number of levels in the discrete variable ",sQuote("y")," is larger than the number of levels in continuous variable ",  sQuote("x")," indicating a possible transposition of the arguments."))
    }
    if(is.factor(x)) {
      stop(paste0("The argument ", sQuote("X"), " is a factor but must be continuous in the Polyserial correlation."))
    }
    if(fast){
      value <- polysFast(x, y, weights, ML=ML)
    }
    else {
     value <- polysSlow(x, y, weights, ML=ML) 
    }
    foundMethod <- TRUE
  }

  if (method == "polychoric") {
    if (fast) {
      value <- polycFast(x, y, w=weights, ML=ML)
    }
    else {
      value <- polycSlow(x, y, w=weights, ML=ML)
    }
    foundMethod <- TRUE
  }

  if (method == "pearson" | method == "spearman") {
    if(fast){
      value <- contCorrFast(x, y, w=weights, method=method)
    }
    else {
      value <- contCorr(x, y, w=weights, method=method)
    }
    foundMethod <- TRUE
  }

  if(!foundMethod) {
    stop(paste0("Could not find method ",sQuote(method0), " see help for available methods."))
  }
  value
}
