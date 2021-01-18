#' convert density function to quantile and quantile density function
#' @export
#' @param quantileCurves n-by-m matrix of quantile curves
#' @param t_vec length m vector contains the common support grid of the quantile curves
quan2den_qd <- function(quantileCurves, t_vec) {

        n = nrow(quantileCurves)
        m = length(t_vec)
        res = list()

        res_list_dSup = sapply(1:n, function(i) {

                quantile_raw = quantileCurves[i, ]
                spline_fit = splinefun(x = t_vec, y = quantile_raw, method = "hyman")
                fitted_quantile = spline_fit(t_vec, deriv = 0)
                fitted_dq = spline_fit(t_vec, deriv = 1)
                fitted_dq_prime = spline_fit(t_vec, deriv = 2)

                ## option1: compute pdf using cubic spline derivative
                spline_fit_cdf = splinefun(x = fitted_quantile, y = t_vec, method = "hyman")
                fitted_pdf = spline_fit_cdf(fitted_quantile, deriv = 1)

                ### correct the fitted_dq = 0 cases,  1/fitted_dq = pdf###
                while (any(fitted_dq==0))  {

                        index = which(fitted_dq == 0)

                        index_plus1 = index + 1
                        index_plus1[index_plus1>m] = m

                        index_miuns1 = index - 1
                        index_miuns1[index_miuns1<1] = 1


                        fitted_dq[index] = fitted_dq[index_plus1]/2 + fitted_dq[index_miuns1]/2

                }
                ## option2: compute pdf using 1/dq
                #fitted_pdf = 1/fitted_dq
                return(cbind(quantile = fitted_quantile, dq = fitted_dq, dq_prime = fitted_dq_prime, fitted_pdf = fitted_pdf))

        }) # the first m rows are quantile functions, middle m rows are quantile density-dq, last m rows are dq'
        res$Qobs = t(res_list_dSup[1:m, ])
        res$qobs = t(res_list_dSup[(m+1):(2*m), ])
        res$qobs_prime = t(res_list_dSup[(2*m+1):(3*m), ])
        res$fobs = t(res_list_dSup[(3*m+1):(4*m), ]) # support is the corresponding value in Qobs

        return(res)

}
