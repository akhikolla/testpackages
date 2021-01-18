#' Perform Frechet Regression with the Wasserstein Distance
#' @export
#' @import stats
#' @import utils
#' @import methods
#' @param rightside_formula a right-side formula
#' @param Xfit_df n-by-p matrix (or dataframe) of predictor values for fitting (do not include a column for the intercept)
#' @param Ytype 'quantile' or 'density'
#' @param Ymat one of the following matrices:
#' \itemize{
#' \item{if Ytype = 'quantile'} Ymat is an n-by-m matrix of the observed quantile functions. Ymat[i, :] is a 1-by-m vector of quantile function values on grid \code{Sup}.
#' \item{if Ytype = 'density'} Ymat is an n-by-m matrix of the observed density functions. Ymat[i, :] is a 1-by-m vector of density function values on grid \code{Sup}.
#' }
#' @param Sup one of the following vectors:
#' \itemize{
#' \item{if Ytype = 'quantile'} Sup is a length m vector - common grid for all quantile functions in Ymat (default: seq(0, 1, length.out = ncol(Ymat))).
#' \item{if Ytype = 'density'} Sup is a length m vector - common grid for all density functions in Ymat (default: seq(0, 1, length.out = ncol(Ymat))).
#' }
#' @return a list containing the following objects:
#' \item{call}{function call}
#' \item{rformula}{\code{rightside_formula}}
#' \item{predictor_names}{names of predictors as the colnames given in the xfit matrix or dataframe.}
#' \item{Qfit}{n-by-m matrix of fitted quantile functions.}
#' \item{xfit}{design matrix in quantile fitting.}
#' \item{Xfit_df}{n-by-p matrix (or dataframe) of predictor values for fitting}
#' \item{Yobs}{a list containing the following matrices:
#' \itemize{
#' \item{Qobs:} {n-by-m matrix of the observed quantile functions.}
#' \item{qobs:} {n-by-m matrix of the observed quantile density functions.}
#' \item{qobs_prime:} {n-by-m matrix of the first derivative of the observed quantile density functions.}
#' \item{fobs:} {n-by-m matrix of the observed density functions.}
#' }}
#' \item{t_vec}{a length m vector - common grid for all quantile functions in Qobs.}
#' @examples
#' data(strokeCTdensity)
#' predictor = strokeCTdensity$predictors
#' dSup = strokeCTdensity$densitySupport
#' densityCurves = strokeCTdensity$densityCurve
#'
#' res1 = wass_regress(rightside_formula = ~., Xfit_df = predictor,
#'  Ytype = 'density', Ymat = densityCurves, Sup = dSup)
#' res2 = wass_regress(rightside_formula = ~ log_b_vol * weight, Xfit_df = predictor,
#'  Ytype = 'density', Ymat = densityCurves, Sup = dSup)
#' @references
#'   \cite{Wasserstein F-tests and confidence bands for the Frechet regression of density response curves, Alexander Petersen, Xi Liu and Afshin A. Divani, 2019}
wass_regress <- function(rightside_formula, Xfit_df, Ytype, Ymat, Sup = NULL) {

        ### ======================  input check  ====================== ###

        if (nargs() < 4) {
                stop("Not enough arguments: rightside_formula, Xfit_df, Ymat and Ytype are required")
        }
        if (is.null(Sup)) Sup = seq(0, 1, length.out = ncol(Ymat))
        if (ncol(Ymat) != length(Sup)) {
                stop("ncol(Ymat) and length(Sup) do not match")
        }
        if (!Ytype %in% c('density', 'quantile')) {
                stop("Ytype must be \"density\" or \"quantile\" ")
        }
        if (Ytype != 'density' & (min(Sup) != 0 | max(Sup) != 1)) {
                stop("Input Sup should be an increasing grid beginning at 0 and ending at 1")
        }
        if (Ytype != 'density' & !all(diff(t(Ytype)) >= 0)) {
                stop("Each row of Ymat should be nondecreasing")
        }
        if (length(Sup) < 25) {
                stop("please give densely observed density/quantile curves, with length(Sup) >= 25 ")
        }

        # Added formula input check
        if (!inherits(rightside_formula, 'formula')) {
                stop("formula argument must be of class \"formula\" ")
        }

        ### ====================== begin function  ====================== ###
        cl <- match.call()

        ### convert Ymat
        if (Ytype == 'density') {
                t_vec = unique(c(seq(0, 0.05, 0.001), seq(0.05, 0.95, 0.05), seq(0.95, 1, 0.001)))
                Yobs = den2Q_qd(densityCurves = Ymat, dSup = Sup, t_vec = t_vec)
        } else {
                Yobs = quan2den_qd(quantileCurves = Ymat, t_vec = Sup)
                t_vec = Sup
        }

        Qmat = Yobs$Qobs

        ### OLS fit of Qmat
        twoside_formula <- as.formula(paste('Qmat', deparse(rightside_formula)))

        # coefficient matrix is a (1+p) x m dimension
        Qmat_lm <- stats::lm(twoside_formula, data = Xfit_df)
        xfit <- model.matrix(Qmat_lm)[ ,-1]
        Qmat_fitted = fitted(Qmat_lm)

        ### quadratic program to make qfitted > 0
        Qmat_fitted = quadraticQ(Qmat_fitted, t_vec)

        ### ======================   return list   ====================== ###

        res = structure(list(
                   call = cl,
                   rformula = rightside_formula,
                   predictor_names = colnames(xfit),

                   Qfit = Qmat_fitted,
                   xfit = xfit,
                   Xfit_df = Xfit_df,
                   Yobs = Yobs,
                   t_vec = t_vec),
                   class = 'WRI')
}
