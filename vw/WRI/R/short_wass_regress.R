#' quick wasserstein regression fitting function with smoothed input
#' @name globalFstat
#' @keywords internal
short_wass_regress <- function(Xfit_df, smoothY){

        Qobs = smoothY$Qobs
        Qmat_lm <- stats::lm(Qobs ~ Xfit_df)
        Qfitted = fitted(Qmat_lm)

        res = structure(list(
                call = NULL,
                lm_res = NULL,
                predictor_names = NULL,

                qfit = NULL,
                ffit = NULL,
                Qfit = Qfitted,

                xfit = Xfit_df,
                Ysmooth = smoothY,
                t_equal = NULL),
                class = 'WRI')
}
