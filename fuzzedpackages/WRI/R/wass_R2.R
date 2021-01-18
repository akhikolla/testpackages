#' Compute Wasserstein Coefficient of Determination
#' @export
#' @references
#' \cite{Frechet regression for random objects with Euclidean predictors, Alexander Petersen and Hans-Georg Müller, 2019}
#' @param wass_regress_res an object returned by the \code{wass_regress} function
#' @return Wasserstein \eqn{R^2}, the Wasserstein coefficient of determination
#' @examples
#' data(strokeCTdensity)
#' predictor = strokeCTdensity$predictors
#' dSup = strokeCTdensity$densitySupport
#' densityCurves = strokeCTdensity$densityCurve
#'
#' res = wass_regress(rightside_formula = ~., Xfit_df = predictor,
#' Ymat = densityCurves, Ytype = 'density', Sup = dSup)
#' wass_r2 = wass_R2(res)
wass_R2 <- function(wass_regress_res) {
        if (class(wass_regress_res) != 'WRI') {
                stop("the first argument should be an object returned by function wass_regress.")
        }
        Qobs = wass_regress_res$Yobs$Qobs
        Qfit = wass_regress_res$Qfit
        t_vec = wass_regress_res$t_vec

        n = nrow(Qobs)
        Qmean = colMeans(Qobs)
        sstotal = sum(sapply(1:n, function(i) fdapace::trapzRcpp(t_vec, (Qobs[i, ] - Qmean)^2)))
        sse = sum(sapply(1:n, function(i) fdapace::trapzRcpp(t_vec, (Qobs[i, ] - Qfit[i, ])^2)))
        wass_r2 = 1 - sse/sstotal
}
