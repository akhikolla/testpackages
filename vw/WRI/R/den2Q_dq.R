#' convert density function to quantile and quantile density function
#' @export
#' @param densityCurves n-by-m matrix of density curves
#' @param dSup length m vector contains the common support grid of the density curves
#' @param t_vec common grid for quantile functions
den2Q_qd <- function(densityCurves, dSup, t_vec) {

        n = nrow(densityCurves)
        m = length(t_vec)
        res = list()
        m_dSup = length(dSup)

        if(min(densityCurves) < 0) {
                stop("Please correct negative or zero probability density estimates.")
        }

        if (length(dSup) < 25) {
                stop("please give densely observed density curves, with length(dSup) >= 25 ")
        }

        res_list_dSup = sapply(1:n, function(i) {

                dens = densityCurves[i, ]
                if (abs(fdapace::trapzRcpp(X = dSup, dens) - 1) > 1e-05) {
                        # warning("Density does not integrate to 1 with tolerance of 1e-5 - renormalizing now.")
                        dens = dens/fdapace::trapzRcpp(X = dSup, Y = dens)
                }

                cdf_raw = fdapace::cumtrapzRcpp(X = dSup, Y = dens)
                spline_fit = splinefun(x = dSup, y = cdf_raw, method = "hyman")
                fitted_cdf = spline_fit(dSup, deriv = 0)
                fitted_pdf = spline_fit(dSup, deriv = 1)
                fitted_pdf_prime = spline_fit(dSup, deriv = 2)

                ### avoid pdf is 0
                if(any(fitted_pdf==0)) {
                        index = which(fitted_pdf == 0)
                        fitted_pdf[index] = min(fitted_pdf[-index])
                        fitted_pdf = fitted_pdf/fdapace::trapzRcpp(X = dSup, Y = fitted_pdf)
                        fitted_cdf = fdapace::cumtrapzRcpp(X = dSup, Y = fitted_pdf)
                }

                return(cbind(cdf = fitted_cdf, pdf = fitted_pdf, pdf_prime = fitted_pdf_prime))
        }) # the first m_dSup rows are CDF, middle m_dSup rows are PDF, last m_dSup rows are PDF'

        Qobs_t = t(apply(res_list_dSup[1:m_dSup, ], MARGIN = 2, FUN = function(fitted_cdf) {
                unique_index = !duplicated(fitted_cdf)
                approx(x = unique(fitted_cdf), y = dSup[unique_index], xout = t_vec,
                       rule = c(2, 2))[[2]]

        }))

        qobs_t = t(sapply(1:n, FUN = function(i) {
                fitted_cdf = res_list_dSup[1:m_dSup, i]
                unique_index = !duplicated(fitted_cdf)
                approx(x = unique(fitted_cdf), y = 1/res_list_dSup[(m_dSup+1):(2*m_dSup), i][unique_index], xout = t_vec,
                       rule = c(2, 2))[[2]]
        }))

        qobs_prime_t = t(sapply(1:n, FUN = function(i) {
                fitted_cdf = res_list_dSup[1:m_dSup, i]
                unique_index = !duplicated(fitted_cdf)
                approx(x = unique(fitted_cdf), y = (-1/res_list_dSup[(m_dSup+1):(2*m_dSup), i]^3 * res_list_dSup[(2*m_dSup+1):(3*m_dSup), i])[unique_index], xout = t_vec,
                       rule = c(2, 2))[[2]]
        }))

        res$Qobs = Qobs_t
        res$qobs = qobs_t
        res$qobs_prime = qobs_prime_t
        res$fobs = 1/res$qobs # support is the corresponding value in Qobs

        return(res)

}
