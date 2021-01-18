#' Confidence Bands for Wasserstein Regression
#' @export
#' @importFrom rlang .data
#' @details This function computes intrinsic confidence bands for \code{Xpred_df} if \code{type} = 'quantile' and density bands if \code{type} = 'density', and visualizes the confidence and/or density bands when \code{figure} = TRUE.
#' @param wass_regress_res an object returned by the \code{wass_regress} function
#' @param Xpred_df k-by-p matrix (or dataframe, or named vector) used for prediction. Note that Xpred_df should have the same column names with Xfit_df used in wass_regress_res
#' @param level confidence level
#' @param delta boundary control value in density band computation. Must be a value in the interval (0, 1/2) (default: 0.01)
#' @param type 'density', 'quantile' or 'both'
#' \itemize{
#' \item{'density': density function bands will be returned (and plotted if \code{figure = TRUE}) }
#' \item{'quantile': quantile function and CDF bands will be returned (and plotted if \code{figure = TRUE}) }
#' \item{'both': three kinds of bands, density function, quantile function and CDF bands will be returned (and plotted if \code{figure = TRUE}) }
#' }
#' @param figure logical; if TRUE, return a sampled plot (default: TRUE)
#' @param fig_num the fig_num-th row of \code{Xpred_df} will be used for visualization of confidence bands. If NULL, then \code{fig_num} is randomly chosen (default: NULL)
#' @return a list containing the following lists:
#' \item{den_list:}{
#' \itemize{
#' \item{fpred:} {k-by-m matrix, predicted density function at Xpred_df.}
#' \item{f_ux:} {k-by-m matrix, upper bound of confidence bands of density functions.}
#' \item{f_lx:} {k-by-m matrix, lower bound of confidence bands of density functions.}
#' \item{Qpred:} {k-by-m matrix, f_lx[i, ], f_ux[i, ] and fpred[i, ] evaluated on Qpred[i, ] vector.}
#' }}
#' \item{quan_list:}{
#' \itemize{
#' \item{Qpred:} {k-by-m matrix of predicted quantile functions.}
#' \item{Q_ux:} {k-by-m matrix of upper bound of quantile functions.}
#' \item{Q_lx:} {k-by-m matrix of lower bound of quantile functions.}
#' \item{t_vec:} {a length m vector - common grid for all quantile functions.}}}
#' \item{cdf_list:}{
#' \itemize{
#' \item{fpred:} {k-by-m matrix, predicted density function.}
#' \item{Fpred:} {k-by-m matrix, predicted cumulative distribution functions.}
#' \item{F_ux:} {k-by-m matrix, upper bound of cumulative distribution functions.}
#' \item{F_lx:} {k-by-m matrix, lower bound of cumulative distribution functions.}
#' \item{Fsup:} {k-by-m matrix, fpred[i, ], F_lx[i, ], F_ux[i, ] and Fpred[i, ] evaluated on Fsup[i, ] vector.}}}
#' @examples
#' alpha = 2
#' beta = 1
#' n = 50
#' x1 = runif(n)
#' t_vec = unique(c(seq(0, 0.05, 0.001), seq(0.05, 0.95, 0.05), seq(0.95, 1, 0.001)))
#' set.seed(1)
#' quan_obs = simulate_quantile_curves(x1, alpha, beta, t_vec)

#' Xfit_df = data.frame(x1 = x1)
#' res = wass_regress(rightside_formula = ~., Xfit_df = Xfit_df,
#'                    Ytype = 'quantile', Ymat = quan_obs, Sup = t_vec)
#' confidence_Band = confidenceBands(res, Xpred_df = data.frame(x1 = c(-0.5,0.5)),
#' type = 'both', fig_num = 2)

#' \donttest{
#' data(strokeCTdensity)
#' predictor = strokeCTdensity$predictors
#' dSup = strokeCTdensity$densitySupport
#' densityCurves = strokeCTdensity$densityCurve
#' xpred = predictor[2:3, ]
#'
#' res = wass_regress(rightside_formula = ~., Xfit_df = predictor,
#' Ytype = 'density', Ymat = densityCurves, Sup = dSup)
#' confidence_Band = confidenceBands(res, Xpred_df = xpred, type = 'density', fig_num = 1)
#' }
confidenceBands <- function(wass_regress_res, Xpred_df, level = 0.95, delta = 0.01, type = 'density', figure = TRUE, fig_num = NULL){

        ### ======================  input check  ====================== ###
        if (class(wass_regress_res) != 'WRI') {
                stop("the first argument should be an object returned by function wass_regress.")
        }
        if (nargs() < 2) {
                stop("Not enough arguments: wass_regress_res and Xpred_df must be provided.")
        }
        if (!type %in% c('density', 'quantile', 'both')) {
                stop("type must be 'density', 'quantile' or 'both'")
        }


        if (is.vector(Xpred_df)) {
                Xpred_df = t(Xpred_df)
                k = 1
        } else { k = nrow(Xpred_df)}

        if (!is.null(fig_num) ) {
                if (fig_num > k) stop("fig_num must be an positive integer less than nrow(Xpred_df)")
        }

        Qobs = wass_regress_res$Yobs$Qobs
        t_vec = wass_regress_res$t_vec
        m = length(t_vec)
        n = nrow(Qobs)

        ### OLS fit of Qmat
        twoside_formula <- as.formula(paste('Qobs', deparse(wass_regress_res$rformula)))
        Qmat_lm <- stats::lm(twoside_formula, data = wass_regress_res$Xfit_df)
        Qfit = fitted(Qmat_lm)
        Qpred = predict(Qmat_lm, newdata = as.data.frame(Xpred_df))
        xpred = model.matrix(wass_regress_res$rformula, data = as.data.frame(Xpred_df))[, -1]

        ### quadratic program to make qfitted > 0
        Qfit = quadraticQ(Qfit, wass_regress_res$t_vec)
        Qpred = quadraticQ(Qpred, wass_regress_res$t_vec)

        ### ====================== 1) begin function  ====================== ###

        Xmat = cbind(1, wass_regress_res$xfit)
        Sigma = t(Xmat) %*% Xmat/n
        Sigma_inv = solve(Sigma)

        Q_res = Qobs - Qfit
        xpred = matrix(xpred, nrow = k)

        qobs = wass_regress_res$Yobs$qobs
        twoside_formula <- as.formula(paste('qobs', deparse(wass_regress_res$rformula)))
        qmat_lm <- stats::lm(twoside_formula, data = wass_regress_res$Xfit_df)
        qfit = fitted(qmat_lm)

        ### to grarantee qfit > 0
        qfit = t(apply(qfit, 1, function(qvec) {
                index = qvec < 0
                qvec[index] = min(qvec[!index])
                return(qvec)
        }))
        qpred = predict(qmat_lm, newdata = as.data.frame(Xpred_df))

        ### to grarantee qpred > 0
        qpred = t(apply(qpred, 1, function(qvec) {
                index = qvec < 0
                qvec[index] = min(qvec[!index])
                return(qvec)
        }))

        qpred_full_support = qpred

        ## ===============        2) compute D(s, t)          =============== ##
        if(type %in% c('density', 'both')) {

                index = which(t_vec > delta & t_vec < (1- delta))
                m_reduce = length(index)

                qpred = matrix(qpred[, index], nrow = k)
                fpred = 1/qpred
                Q_res_den = Q_res[ ,index]
                qobs = qobs[, index]
                qfit = qfit[, index]
                q_res = qobs - qfit
                t_reduce = t_vec[index]

                # Rcpp::sourceCpp("src/Rcpp_D_4cell.cpp")
                rcpp_4 = Rcpp_D_4cell(Xmat = Xmat, Q_res = Q_res_den, q_res = q_res)

                qobs_prime = wass_regress_res$Yobs$qobs_prime[, index]
                twoside_formula <- as.formula(paste('qobs_prime', deparse(wass_regress_res$rformula)))
                qmat_prime_lm <- stats::lm(twoside_formula, data = wass_regress_res$Xfit_df)
                q_prime_pred = predict(qmat_prime_lm, newdata = as.data.frame(Xpred_df)) # k x m_reduce

        }

        if (type %in% c('quantile', 'both')) {
                # Rcpp::sourceCpp("src/Rcpp_D_cell.cpp")
                D_cell_rcpp =  Rcpp_D_cell(Xmat = Xmat, Q_res = Q_res[, c(-1, -m)])

        }

        ## =================================================================== ##
        ##                           quantile band                             ##
        ## =================================================================== ##
        ## =========    for l = 1:k  3.A) compute m_alpha and se   =========== ##
        if (type %in% c('quantile', 'both')) {

                m_alpha = rep(0, k)
                se = matrix(0, nrow = k, ncol = m)
                C_x = matrix(0, nrow = m-2, ncol = m-2)

                for (l in 1:k) {
                        x_star = Sigma_inv %*% c(1, xpred[l, ])
                        for (i in 1:(m - 2)) {
                                for (j in 1:i){
                                        C_x[i, j] = t(x_star) %*% D_cell_rcpp[[(i-1)*(m-2) + j]] %*% x_star
                                        C_x[j, i] = C_x[i, j]
                                }
                        }

                        # --------------------------   compute m_alpha    ------------------------- #
                        C_x_diag = sqrt(diag(C_x))
                        Ruu_mat = matrix(rep(C_x_diag, ncol(C_x)), nrow = ncol(C_x))
                        C_x_stand = C_x/Ruu_mat/t(Ruu_mat)

                        GaussianProcess = mvtnorm::rmvnorm(n = 1000, mean = rep(0, m-2) ,sigma = C_x_stand)
                        sequence_max = Rfast::rowMaxs(abs(GaussianProcess), value = TRUE)

                        m_alpha[l] = quantile(sequence_max, probs = level)
                        # -----------   compute se = n^{-1/2} * (x' * hat{C}(t,t) * x)  ----------- #
                        C_x_diag = c(C_x_diag[1], C_x_diag, C_x_diag[m-2])
                        se[l, ] = C_x_diag * sqrt(1/n)
                }
                ## ==================   4.A) compute Q_lx and Q_ux      ================== ##

                Q_lx = matrix(0, nrow = k, ncol = m)
                Q_ux = matrix(0, nrow = k, ncol = m)


                for (i in 1:k){
                        Qi_lx = Qpred[i, ] - m_alpha[i] * se[i, ]
                        Qi_ux = Qpred[i, ] + m_alpha[i] * se[i, ]

                        if (any(diff(Qi_lx) < 0)){
                                Qrefit = CVXR::Variable(m)
                                obj = sum(((Qi_lx - Qrefit))^2)
                                prob = CVXR::Problem(CVXR::Minimize(obj), list(CVXR::diff(Qrefit) >= 0, Qpred[i, ] - Qrefit >= 0))
                                result = CVXR::psolve(prob)
                                Qrefit = result$getValue(Qrefit)
                                Qi_lx = Qrefit
                        }

                        if (any(diff(Qi_ux) < 0)){
                                Qrefit = CVXR::Variable(m)
                                obj = sum(((Qi_ux - Qrefit))^2)
                                prob = CVXR::Problem(CVXR::Minimize(obj), list(CVXR::diff(Qrefit) >= 0,  Qrefit - Qpred[i, ] >= 0))
                                result = CVXR::psolve(prob)
                                Qrefit = result$getValue(Qrefit)
                                Qi_ux = Qrefit
                        }
                        Q_lx[i, ] = Qi_lx
                        Q_ux[i, ] = Qi_ux
                }
                quan_list = list(Q_lx = Q_lx,
                                 Qpred = Qpred,
                                 Q_ux = Q_ux,
                                 t_vec = t_vec)

                ## ================          add CDF bands           ================ ##
                cdf_list = list(
                        fpred = 1/qpred_full_support,
                        Fpred = matrix(rep(t_vec, k), nrow = k, byrow = TRUE),
                        Fsup = Qpred,
                        F_lx = t(sapply(1:k, function(i) {approx(x = Q_ux[i, ], y = t_vec, xout = Qpred[i, ], rule = 2)$y})),
                        F_ux = t(sapply(1:k, function(i) {approx(x = Q_lx[i, ], y = t_vec, xout = Qpred[i, ], rule = 2)$y})))
        } else {
                quan_list = NULL
                cdf_list = NULL
        }


        ## =================================================================== ##
        ##                           density band                              ##
        ## =================================================================== ##
        if (type %in% c('density', 'both')) {
                ### ============       for l = 1:k  3.B1) compute tau(s, t)     ========== ###

                l_alpha = rep(0, k)
                se = matrix(0, nrow = k, ncol = m_reduce)
                Tau_x = matrix(0, m_reduce, m_reduce)

                for (l in 1:k) {
                        ct_2_m = rbind(q_prime_pred[l, ], qpred[l, ]) %*% diag(1/qpred[l, ]^3)
                        x_star = Sigma_inv %*% c(1, xpred[l, ])

                        for (i in 1:m_reduce){
                                for (j in 1:i){
                                        Tau_x[i, j] = t(ct_2_m[, i]) %*% matrix(c(
                                                t(x_star) %*% rcpp_4$D_00[[(i-1)*m_reduce + j]] %*% x_star,
                                                t(x_star) %*% rcpp_4$D_10[[(i-1)*m_reduce + j]] %*% x_star,
                                                t(x_star) %*% rcpp_4$D_01[[(i-1)*m_reduce + j]] %*% x_star,
                                                t(x_star) %*% rcpp_4$D_11[[(i-1)*m_reduce + j]] %*% x_star
                                        ), nrow = 2, ncol = 2) %*% ct_2_m[, j]
                                        Tau_x[j, i] =  Tau_x[i, j]
                                }
                        }
                        ### ============       for l = 1:k  3.B2) compute l_alpha       ============ ###

                        ### from covariance Tau_x to correlation Tau_x_c
                        Tau_x_diag = sqrt(diag(Tau_x))
                        Ruu_mat = matrix(rep(Tau_x_diag, m_reduce), m_reduce, m_reduce)
                        Rvv_mat = t(Ruu_mat)
                        Tau_x_c = Tau_x/(Ruu_mat*Rvv_mat) # standardize C_x have 1 on diagonal

                        GaussianProcess = mvtnorm::rmvnorm(n = 1000, sigma = Tau_x_c)
                        sequence_max = Rfast::rowMaxs(abs(GaussianProcess), value = TRUE)

                        l_alpha[l] = quantile(sequence_max, probs = level)
                        # -----------   compute se = n^{-1/2} * (x' * hat{C}(t,t) * x)  ----------- #
                        se[l, ] = Tau_x_diag * sqrt(1/n)
                }

                ### ========    4.B)  compute f_lx and f_ux and form res_list.   ========= ###

                if (k == 1) {
                        f_lx = fpred - l_alpha * se
                        f_ux = fpred + l_alpha * se
                } else {
                        f_lx = fpred - diag(l_alpha) %*% se
                        f_ux = fpred + diag(l_alpha) %*% se
                }

                f_lx[f_lx<0] = 0

                den_list = list(
                        fpred = matrix(fpred, nrow = k),
                        f_lx = f_lx,
                        f_ux = f_ux,
                        Qpred = matrix(Qpred[, index], nrow = k)
                )
        } else {
                den_list = NULL
        }

        ## ==================      5) plot density bands      ================== ##
        if (figure) {
                if (is.null(fig_num)) {
                        fig_num = sample(1:k, size = 1)
                }

                if (type %in% c('quantile', 'both')) {
                        plotdf = data.frame(lower_bound = Q_lx[fig_num, ],
                                            quantile = Qpred[fig_num, ],
                                            upper_bound = Q_ux[fig_num, ],
                                            t = t_vec)

                        p1 =  ggplot2::ggplot(data = plotdf, ggplot2::aes(x = t)) +
                                ggplot2::geom_line(ggplot2::aes(y = .data$quantile), col = 'slateblue1', size = 0.8) +
                                ggplot2::geom_ribbon(ggplot2::aes(ymin = .data$lower_bound, ymax = .data$upper_bound), fill = 'grey70', alpha = 0.8)+
                                ggplot2::ylab('Quantile Function')+
                                ggplot2::theme_linedraw()

                        plotdf = data.frame(lower_bound = cdf_list$F_lx[fig_num, ],
                                            Fsup = cdf_list$Fsup[fig_num, ],
                                            upper_bound = cdf_list$F_ux[fig_num, ],
                                            Fpred = cdf_list$Fpred[fig_num, ])

                        p3 = ggplot2::ggplot(data = plotdf, ggplot2::aes(x = .data$Fsup)) +
                                ggplot2::geom_line(ggplot2::aes(y = .data$Fpred), col = 'slateblue1', size = 0.8) +
                                ggplot2::geom_ribbon(ggplot2::aes(ymin = .data$lower_bound, ymax = .data$upper_bound), fill = 'grey70', alpha = 0.8)+
                                ggplot2::ylab('CDF') +
                                ggplot2::xlab('Support')+
                                ggplot2::theme_linedraw()

                }
                if (type %in% c('density', 'both')) {
                        plotdf = data.frame(lower_bound = f_lx[fig_num, ],
                                            density = den_list$fpred[fig_num, ],
                                            upper_bound = f_ux[fig_num, ],
                                            support = den_list$Qpred[fig_num, ])
                        plotdf_pdf_full = data.frame(support = Qpred[fig_num, ],
                                                     density = 1/qpred_full_support[fig_num, ])

                        p2 =  ggplot2::ggplot(data = plotdf, ggplot2::aes(x = .data$support)) +
                                ggplot2::geom_ribbon(ggplot2::aes(ymin = .data$lower_bound, ymax = .data$upper_bound), fill = 'grey70', alpha = 0.8)+
                                ggplot2::geom_line(data = plotdf_pdf_full, ggplot2::aes(x = .data$support, y = .data$density), col = 'slateblue1', size = 0.8) +
                                ggplot2::ylab('PDF') +
                                ggplot2::xlab('Support') +
                                ggplot2::theme_linedraw()


                }
                if (type == 'density') {
                        show(p2 + ggplot2::ggtitle(paste0('Confidence Band of Density Function for Object ', fig_num)))
                } else if (type == 'quantile') {
                        pp = gridExtra::grid.arrange(p1, p3, top = paste0('Confidence Bands for Object ', fig_num), nrow = 1)
                } else if (type == 'both') {
                        pp = gridExtra::grid.arrange(p1, p3, p2, top = paste0('Confidence Bands for Object ', fig_num), nrow = 1)
                }
        }

        return (list(
                quan_list = quan_list,
                cdf_list = cdf_list,
                den_list = den_list
        ))

}








