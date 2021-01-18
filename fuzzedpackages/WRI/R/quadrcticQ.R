#' An internal function to do quadratic program in order to make Qfitted nondescreasing
#' @name quadraticQ
#' @keywords internal
quadraticQ <- function(Qmat, t) {

        Qmat_diff = t(diff(t(Qmat)))

        t_diff = diff(t)
        t_diff_m = 0.5*(c(0, t_diff) + c(t_diff, 0))

        if (min(Qmat_diff) < 0) {

                message('use quadratic programming to make fitted quantile density >= 0')

                refit_index = which(Rfast::rowMins(Qmat_diff, value = TRUE) < 0)
                for (i in refit_index) {
                        ### solve constrained optimization
                        m = length(Qmat[i, ])
                        Qrefit = CVXR::Variable(m)
                        obj = sum(((Qmat[i, ] - Qrefit)*t_diff_m)^2)
                        prob = CVXR::Problem(CVXR::Minimize(obj), list(CVXR::diff(Qrefit) >= 0))
                        result = CVXR::psolve(prob)
                        Qrefit = result$getValue(Qrefit)

                        Qmat[i, ] = Qrefit
                }
        }

        return(Qmat)

}
