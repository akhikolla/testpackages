#' Summary Function of Wasserstein Regression Model
#' @export
#' @param object an object returned by the \code{wass_regress} function
#' @param ... further arguments passed to or from other methods.
#' @return a list containing the following fields:
#' \item{call}{function call of the Wasserstein regression}
#' \item{r.square}{Wasserstein \eqn{R^2}, the Wasserstein coefficient of determination}
#' \item{global_wasserstein_F_stat}{Wasserstein global F test statistic from the Satterthwaite method}
#' \item{global_F_pvalue}{p value of global F test}
#' \item{global_wasserstein_F_df}{degrees of freedom of satterthwaite approximated sampling distribution used in global F test}
#' \item{partial_F_table}{Partial F test for individual effects}
#' @examples
#' data(strokeCTdensity)
#' predictor = strokeCTdensity$predictors
#' dSup = strokeCTdensity$densitySupport
#' densityCurves = strokeCTdensity$densityCurve
#'
#' res <- wass_regress(rightside_formula = ~., Xfit_df = predictor,
#' Ymat = densityCurves, Ytype = 'density', Sup = dSup)
#' summary(res)
summary.WRI <- function(object, ...) {
        digits = 3
        ans = list()
        ans$call = object$call
        ### ================== compute Wasserstain R^2 ===================== ###
        ans$r.square = round(wass_R2(object), digits = digits)

        ### ==================   compute global F test ===================== ###
        globalFres = globalFtest(object, alpha = 0.05, permutation = FALSE, numPermu = 200, bootstrap = FALSE, numBoot = 200)
        ans$global_F_pvalue = globalFres$summary_df$p_value[2]
        ans$global_wasserstein_F_stat = round(globalFres$wasserstein.F_stat, digits = digits)
        ans$global_wasserstein_F_df =  round(globalFres$chisq_df, digits = digits)

        ### ==================   compute partial F test ==================== ###
        p = ncol(object$xfit)
        coef_table = sapply(1:p, function(i) {
                reduced_model = short_wass_regress(Xfit_df = object$xfit[, -i], smoothY = object$Yobs)
                partialFres = partialFtest(object, reduced_res = reduced_model, alpha = 0.05)
                return(c(partialFres$statistic[1], partialFres$p_value))
        })
        coef_table = t(coef_table)
        rownames(coef_table) = object$predictor_names
        colnames(coef_table) = c('F-stat ', ' p-value(truncated) ', ' p-value(satterthwaite)')
        ans$partial_F_table = round(coef_table, digits = digits)

        class(ans) = 'summary.WRI'
        return(ans)
}
