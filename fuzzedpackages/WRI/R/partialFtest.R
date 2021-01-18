#' partial F test for Wasserstein regression
#' @export
#' @details two methods used to compute p value using asymptotic distribution of F statistic
#' \itemize{
#' \item{truncated:}{ asymptotic inference, p-value is obtained by truncating the infinite summation of eigenvalues into the first K terms, where the first K terms explain more than 99.99\% of the variance.}
#' \item{satterthwaite:}{ asymptotic inference, p-value is computed using Satterthwaite approximation method of mixtures of chi-square.}}
#' @param reduced_res a reduced model list returned by the \code{wass_regress} function
#' @param full_res a full model list returned by the \code{wass_regress} function
#' @param alpha type one error rate
#' @return a dataframe containing the following columns:
#' \item{method}{methods used to compute p value, see details}
#' \item{statistic}{the test statistics}
#' \item{critical_value}{critical value}
#' \item{p_value}{p value of global F test}
#' @examples
#' data(strokeCTdensity)
#' predictor = strokeCTdensity$predictors
#' dSup = strokeCTdensity$densitySupport
#' densityCurves = strokeCTdensity$densityCurve
#'
#' full_res <- wass_regress(rightside_formula = ~., Xfit_df = predictor,
#'  Ymat = densityCurves, Ytype = 'density', Sup = dSup)
#' reduced_res <- wass_regress(~ log_b_vol + b_shapInd + midline_shift + B_TimeCT, Xfit_df = predictor,
#'  Ymat = densityCurves, Ytype = 'density', Sup = dSup)
#' partialFtable = partialFtest(reduced_res, full_res, alpha = 0.05)

partialFtest <- function(reduced_res, full_res, alpha = 0.05) {
        if (class(full_res) != 'WRI') {
                stop("the first argument should be an object returned by function wass_regress.")
        }
        if (class(reduced_res) != 'WRI') {
                stop('the second argument should be a list returned by function wass_regress')
        }
        p = ncol(full_res$xfit)
        q = ncol(reduced_res$xfit)
        m = length(full_res$t)
        n = nrow(full_res$xfit)
        t_vec = full_res$t_vec
        diff_t = diff(t_vec)
        Qobs = full_res$Yobs$Qobs
        Qfit_full = full_res$Qfit

        ### ======================  F test statistic  ====================== ###
        Qmean = colMeans(Qobs)
        F_full = sum(sapply(1:n, function(i) fdapace::trapzRcpp(t_vec, (Qmean - full_res$Qfit[i, ])^2)))
        F_reduced = sum(sapply(1:n, function(i) fdapace::trapzRcpp(t_vec, (Qmean - reduced_res$Qfit[i, ])^2)))
        partialFstat = F_full - F_reduced

        ### ==============   Asymptotic r dimensional kernel   ==============###
        r = p - q
        combinedX = cbind(reduced_res$xfit, full_res$xfit)
        X = combinedX[, !duplicated(combinedX, MARGIN = 2)]
        Xc = scale(X, center = TRUE, scale = FALSE)

        ### -------------    reduce m to speed up computation  ------------- ###
        Q_residual = Qobs - Qfit_full
        sigmaXX = cov(X)
        sigmaZY = matrix(sigmaXX[(q+1):p, 1:q], nrow = r, ncol = q)
        # sigmaZY = sigmaXX[(q+1):p, 1:q]
        sigmaYY = sigmaXX[1:q, 1:q]
        sigmaZZ = sigmaXX[(q+1):p, (q+1):p]

        J = cbind(-sigmaZY %*% solve(sigmaYY), diag(rep(1, r)))
        sigmaZ_Y = sigmaZZ - sigmaZY %*% solve(sigmaYY)  %*% t(sigmaZY)
        left_mat = solve(expm::sqrtm(sigmaZ_Y)) %*% J

        # Rcpp::sourceCpp("R/kernel_partialF.cpp")
        r_dim_kernal = kernel_partialF(Xc = Xc, Q_res = Q_residual, left_mat = left_mat, r = r)


        trace_r_dim_kernal = sum(diag(r_dim_kernal)) *  1/(m-1)
        HS_r_dim_kernal = sum(diag(t(r_dim_kernal) %*% r_dim_kernal)) * 1/(m-1) * 1/(m-1)

        ### ======================   truncated method ====================== ###
        eigenvalues = eigen(r_dim_kernal)$values/(m - 1)
        eigenvalues = eigenvalues[eigenvalues > 0]
        trun_index = cumsum(eigenvalues)/sum(eigenvalues) < 0.999
        Keigenvalue = eigenvalues[trun_index]
        repN = 500
        chisquare = matrix(rchisq(repN * length(Keigenvalue), df = r), nrow = repN, ncol = length(Keigenvalue))
        empiricalF = chisquare %*% Keigenvalue
        critical_value_trun = quantile(empiricalF, probs = 1 - alpha)
        pValue_trun = (sum(empiricalF > partialFstat) + 1)/(repN + 1)

        ### =================   satterthwaite method   ===================== ###

        a = HS_r_dim_kernal/trace_r_dim_kernal
        m_chisq = r * trace_r_dim_kernal / a
        criticalValue_satt = qchisq(1 - alpha, m_chisq) * a;
        pValue_satt = pchisq(partialFstat/a, m_chisq, lower.tail = F)

        ### ========================   return list   ======================= ###
        res_df = data.frame(method = c('truncated', 'satterthwaite') , statistic = c(partialFstat, partialFstat),
                            critical_value = c(critical_value_trun, criticalValue_satt),
                            p_value = c(pValue_trun, pValue_satt))
}
