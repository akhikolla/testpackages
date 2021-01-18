#' Simulate quantile curves
#' @export
#' @description This function simulates quantile curves used as a toy example
#' @param x1 n-by-1 predictor vector
#' @param alpha parameter in location transformation
#' @param beta parameter in variance transformation
#' @param t_vec a length m vector - common grid for all quantile functions
#' @return {quan_obs} {n-by-m matrix of quantile functions}
#' @examples
#' alpha = 2
#' beta = 1
#' n = 100
#' x1 = runif(n)
#' t_vec = unique(c(seq(0, 0.05, 0.001), seq(0.05, 0.95, 0.05), seq(0.95, 1, 0.001)))
#' quan_obs = simulate_quantile_curves(x1, alpha, beta, t_vec)
#' @references
#'   \cite{Wasserstein F-tests and confidence bands for the Frechet regression of density response curves, Alexander Petersen, Xi Liu and Afshin A. Divani, 2019}
simulate_quantile_curves <- function(x1, alpha, beta, t_vec) {

        n = length(x1)
        m = length(t_vec)
        mu0 = 0
        sigma0 = 2
        p1 = pnorm(2.5)
        p0 = pnorm(-2.5)

        ### ====================== 1) marginal means  ====================== ###

        Qstar = qnorm((p1-p0)*t_vec + p0)
        qstar = (p1 - p0)/dnorm(Qstar)
        qstar_prime = qstar^3 * Qstar * dnorm(Qstar)/(p1-p0)
        fstar = dnorm(Qstar)/(p1-p0)

        ### =================== 2) conditional means  ====================== ###

        ### - 2.1) intercept and slope used in marginal to conditional mean - ###
        mu_x_vec = mu0 + alpha*x1
        sigma_x_vec = sigma0 + beta*x1 #n-by-1

        ### -------------- 2.2) marginal to conditional mean --------------- ###
        Qmean = matrix(rep(mu_x_vec, m), nrow = n) + diag(sigma_x_vec)%*%matrix(rep(Qstar, n), nrow = n, ncol = m, byrow = TRUE)

        a = runif(n)*4 - 2
        b = runif(n)/2 + 0.75
        Q_obs = matrix(rep(a, m), nrow = n) + diag(b)%*%Qmean
        return(Q_obs)
        }

