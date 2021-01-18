#' global F test for Wasserstein regression
#' @export
#' @details four methods used to compute p value of global F test
#' \itemize{
#' \item{truncated:}{ asymptotic inference,  p-value is obtained by truncating the infinite summation of eigenvalues into the first K terms, where the first K terms explain more than 99.99\% of the variance.}
#' \item{satterthwaite:}{ asymptotic inference, p-value is computed using Satterthwaite's approximation method of mixtures of chi-square.}
#' \item{permutation:}{ resampling technique; Wasserstein SSR is used as the F statistic.}
#' \item{bootstrap:}{ resampling technique; Wasserstein SSR is used as the F statistic.}}
#' @param wass_regress_res an object returned by the \code{wass_regress} function
#' @param alpha type one error rate
#' @param permutation logical; perform permutation global F test (default: FALSE)
#' @param numPermu number of permutation samples if permutation = TRUE
#' @param bootstrap logical; bootstrap global F test (default: FALSE)
#' @param numBoot number of bootstrap samples if bootstrap = TRUE
#' @return a list containing the following fields:
#' \item{wasserstein.F_stat}{the Wasserstein F statistic value in Satterthwaite method .}
#' \item{chisq_df}{the degree of freedom of the null chi-square distribution.}
#' \item{summary_df}{a dataframe containing the following columns:}
#' \itemize{
#' \item{method:} {methods used to compute p value, see details}
#' \item{statistic:} {the test statistics}
#' \item{critical_value:} {critical value}
#' \item{p_value:} {p value of global F test}}
#' @examples
#' data(strokeCTdensity)
#' predictor = strokeCTdensity$predictors
#' dSup = strokeCTdensity$densitySupport
#' densityCurves = strokeCTdensity$densityCurve
#'
#' res = wass_regress(rightside_formula = ~., Xfit_df = predictor,
#'  Ytype = 'density', Ymat = densityCurves, Sup = dSup)
#' globalF_res = globalFtest(res, alpha = 0.05, permutation = TRUE, numPermu = 200)

globalFtest <- function(wass_regress_res, alpha = 0.05, permutation = FALSE, numPermu = 200, bootstrap = FALSE, numBoot = 200) {

        Qobs = wass_regress_res$Yobs$Qobs
        Qfit = wass_regress_res$Qfit
        n = nrow(Qobs)
        t_vec = wass_regress_res$t_vec
        m = length(t_vec)
        diff_t = diff(t_vec)

        ### ======================  F test statistic  ====================== ###
        Qmean = colMeans(Qobs)
        Fstat = sum(sapply(1:n, function(i) fdapace::trapzRcpp(t_vec, (Qmean - Qfit[i, ])^2)))

        ### covariance operator C_Q(s, t)
        resi_mat = Qobs - Qfit
        C_Q = cov(resi_mat)
        trace_C_Q =  fdapace::trapzRcpp(t_vec, diag(C_Q))
        delta_sdelta_t = diff_t %*% t(diff_t)
        HSnorm_C_Q = sum(C_Q[-1, -1]^2*delta_sdelta_t)

        ### ======================   truncated method ====================== ###

        eigenvalues = eigen(C_Q)$values/(m - 1)
        eigenvalues = eigenvalues[eigenvalues > 0]
        trun_index = cumsum(eigenvalues)/sum(eigenvalues) < 0.9999
        Keigenvalue = eigenvalues[trun_index]
        repN = 500
        chisquare = matrix(rchisq(repN * length(Keigenvalue), df = ncol(wass_regress_res$xfit)), nrow = repN, ncol = length(Keigenvalue))
        empiricalF = chisquare %*% Keigenvalue
        critical_value_trun = quantile(empiricalF, probs = 1 - alpha)
        pValue_trun = (sum(empiricalF > Fstat) + 1)/(repN + 1)

        ### =================   satterthwaite method   ===================== ###
        a = HSnorm_C_Q/trace_C_Q
        m_chisq = ncol(wass_regress_res$xfit) * trace_C_Q / a
        criticalValue_satt = qchisq(1 - alpha, m_chisq) * a;
        pValue_satt = pchisq(Fstat/a, m_chisq, lower.tail = F)

        ### ========================   return list   ======================= ###
        res_df = data.frame(method = c('truncated', 'satterthwaite') , statistic = c(Fstat, Fstat),
                            critical_value = c(critical_value_trun, criticalValue_satt),
                            p_value = c(pValue_trun, pValue_satt))

        ### ====================   permutation method   ==================== ###
        if (permutation) {
                Fstat_vec = sapply(1:numPermu, function(i) {
                        set.seed(i)
                        index = sample(1:n, size = n, replace = FALSE)
                        fstat = globalFstat(wass_regress_res$xfit, Qobs[index, ], t_vec)
                        return(fstat)
                })
                pvalue_permu = (sum(Fstat_vec > Fstat) + 1)/(numPermu + 1)
                critical_value_permu = quantile(Fstat_vec, probs = 1 - alpha)
                res_df = rbind(res_df, data.frame(
                        method = 'permutation' , statistic = Fstat,
                        critical_value = critical_value_permu,
                        p_value = pvalue_permu
                ))
        }

        ### ====================   bootstrap method   ==================== ###
        if (bootstrap) {
                Fstat_vec = sapply(1:numBoot, function(i) {
                        set.seed(i)
                        index = sample(1:n, size = n, replace = TRUE)
                        fstat = globalFstat(wass_regress_res$xfit, Qobs[index, ], t_vec)
                        return(fstat)
                })
                pvalue_boot = (sum(Fstat_vec > Fstat) + 1)/(numBoot + 1)
                critical_value_boot = quantile(Fstat_vec, probs = 1 - alpha)
                res_df = rbind(res_df, data.frame(
                        method = 'bootstrap' , statistic = Fstat,
                        critical_value = critical_value_boot,
                        p_value = pvalue_boot
                ))
        }

        rownames(res_df) = NULL
        res_list = list(summary_df = res_df,
                        wasserstein.F_stat = Fstat/a,
                        chisq_df = m_chisq)
        return(res_list)
}

