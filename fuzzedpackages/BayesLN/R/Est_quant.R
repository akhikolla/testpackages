#' Integral of the target functional to minimize
#'
#' Function that implements the integrand of the target functional to minimize in order to
#' obtain hyperparameters that produces optimal estimates in the frequentist context.
#'
#' It is implemented exploiting the infinite sum representation of the estimator.
#'
#'
#'@keywords internal



g_V <- function(Y, lambda, delta, gamma, quant, sigma2, n, nreg, rel_tol = 1e-4) {
  l <- lambda - (n / 2 - nreg / 2)
  w <- nreg / n
  ga <- gamma / sqrt(w)
  beta <- qnorm(p = quant, mean = 0, sd = 1) / sqrt(w)
  d <- sqrt(Y + delta ^ 2) * sqrt(w)
  MGF <-  SMNG_MGF(r = 1, mu = 0, lambda = l, delta = d,
                  gamma = ga, beta = beta, inf_sum = TRUE, rel_tol = rel_tol)
  target <- (MGF - exp(qnorm(quant) * sqrt(sigma2) - (3 * sigma2 * w) / (2))) ^ 2 *
    dgamma(x = Y, shape = (n - nreg) / 2, rate = 1 / (2 * sigma2))
  target <- ifelse(is.finite(target), yes = target, no = 0)
  return(target)
}


#' Vectorization of the function \code{\link{g_V}}
#'
#' It allows to use the integrand of the target functional to
#' minimize as the argument of the function \code{\link{integrate}}.
#'
#'
#'@keywords internal

g_V_vec <- function(Y, lambda, delta, gamma, quant, sigma2, n, nreg=nreg, rel_tol = 1e-4) {
  vec <- sapply(X = Y, function(X) g_V(Y = X, lambda = lambda, delta = delta, gamma = gamma,
                                        quant = quant, sigma2 = sigma2, n = n, nreg=nreg, rel_tol = rel_tol))
  return(vec)
}



#'Target functional to minimize with respect to gamma
#'
#'Function that evaluate the target functional to minimize
#'allowing for a different value of the parameter \code{gamma}.
#'
#'@keywords internal

functional_gamma <- function(x, n, quant, delta, sigma2, lambda, nreg=1, rel_tol = 1e-4) {
  E_gY2 <- integrate(f = g_V_vec, lower = 0, upper = Inf, n = n, lambda = lambda, delta = delta,
                     gamma = x, quant = quant, sigma2 = sigma2,  nreg=nreg, rel.tol = rel_tol)$value
  return(E_gY2)
}



#'Target functional to minimize with respect to delta
#'
#'Function that evaluate the target functional to minimize allowing
#'for a different value of the parameter \code{delta}.
#'
#'@keywords internal


functional_delta <- function(x, n, quant, gamma, sigma2,lambda, nreg=1, rel_tol = 1e-4) {
  E_gY2 <- integrate(f = g_V_vec, lower = 0, upper = Inf, n = n, lambda = lambda, delta = x,
                     gamma = gamma, quant = quant, sigma2 = sigma2,  nreg=nreg, rel.tol = rel_tol)$value
  return(E_gY2)
}




#'Bayesian estimate of the log-normal quantiles
#'
#'This function produces an estimate for the log-normal distribution quantile of fixed level \code{quant}.
#'
#'@param x Vector of data used to estimate the quantile.
#'@param quant Number between 0 and 1 that indicates the quantile of interest.
#'@param method String that indicates the prior setting to adopt. Choosing \code{"weak_inf"}
#'a weakly informative prior setting is adopted, whereas selecting
#' \code{"optimal"} the hyperparameters are fixed trough a numerical optimization algorithm
#' aimed at minimizing the frequentist MSE.
#' @param guess_s2 Specification of a guess for the variance if available. If not, the sample estimate is used.
#' @param x_transf Logical. If \code{TRUE}, the \code{x} vector is assumed already log-transformed.
#' @param CI Logical. With the default choice \code{TRUE}, the posterior credibility interval is computed.
#' @param alpha_CI Level of alpha that determines the credibility (1-\code{alpha_CI}) of the posterior interval.
#' @param type_CI String that indicates the type of interval to compute: \code{"two-sided"} (default),
#'\code{"UCL"} (i.e. Upper Credible Limit) for upper one-sided intervals  or \code{"LCL"} (i.e. Lower
#'Credible Limit) for lower one-sided intervals.
#'@param method_CI String that indicates if the limits should be computed through the logSMNG
#'quantile function \code{\link{qlSMNG}} (option \code{"exact"}, default), or by randomly generating a sample
#'(\code{"simulation"}) using the function \code{\link{rlSMNG}}.
#'@param rel_tol_CI Level of relative tolerance required for the \code{integrate} procedure or for the infinite sum.
#' Default set to \code{1e-5}.
#'@param nrep_CI Number of simulations in case of \code{method="simulation"}.
#'
#'@details
#' The function allows to carry out Bayesian inference for the unconditional quantiles of a sample that is assumed log-normally distributed.
#'
#' A generalized inverse Gaussian prior is assumed for the variance in the log scale \eqn{\sigma^2}, whereas a
#' flat improper prior is assumed for the mean in the log scale \eqn{\xi}.
#'
#' Two alternative hyperparamters setting are implemented (choice controlled by the argument \code{method}): a weakly
#' informative proposal and an optimal one.
#'
#'
#'
#'@return
#' The function returns the prior parameters and their posterior values, summary statistics of the log-scale parameters and the estimate of the specified quantile:
#' the posterior mean and variance are provided by default. Moreover, the user can control the computation of posterior intervals.
#'
#'@source
#'
#'Gardini, A., C. Trivisano, and E. Fabrizi. \emph{Bayesian inference for quantiles of the log-normal distribution.} Biometrical Journal (2020).
#'
#'
#' @examples
#'library(BayesLN)
#'data("EPA09")
#'LN_Quant(x = EPA09, quant = 0.95, method = "optimal", CI = FALSE)
#'LN_Quant(x = EPA09, quant = 0.95, method = "weak_inf",
#'         alpha_CI = 0.05, type_CI = "UCL")
#'
#' @export



LN_Quant <- function(x, quant, method = "weak_inf", x_transf = TRUE, guess_s2 = NULL,CI = TRUE,
                     alpha_CI = 0.05, type_CI = "two-sided", method_CI = "exact",
                     rel_tol_CI = 1e-5, nrep_CI = 1e6){
  if(is.vector(x) == FALSE)
    stop("x must be a vector")
  if(quant > 1 | quant < 0)
    stop("quant must be between 0 and 1")
  if (x_transf == FALSE){
    x <- log(x)
  } else {
    x <- x
  }
  mu <- mean(log(x))
  s2 <- var(log(x))
  n <- length(x)
  beta <- qnorm(p = quant,mean = 0,sd = 1) * sqrt(n)
  if(method == "weak_inf"){
    l <- -(n / 2 - 0.5)
    g <- (3 / sqrt(n)) * sqrt(n)
    d <- sqrt(sum((log(x) - mu) ^ 2) + 0.1 ^ 2) / sqrt(n)
    par <- round(matrix(c(0, 0.01, 3 / sqrt(n),
                          NA, NA, l, d, g, mu, beta), ncol = 5, byrow = T), 3)
  }else if(method=="optimal"){
    l<--(n/2-0.5)
    if(quant==0.5){
      stop("The 'numerical' method cannot be used for the median: use the 'weak_inf' method")
    }else if(quant > 0.5){
      g0 <- (3 / sqrt(n))
      delta <- ifelse(is.numeric(guess_s2), yes = guess_s2, no = 1)
      s2 <- ifelse(is.numeric(guess_s2), yes = guess_s2, no = s2)
      gamma <- optimx::optimx(fn = functional_gamma, par = g0 + 0.7,lower = g0, method = "nlminb", n = n,
                              quant = quant, delta = delta, sigma2 = s2, lambda = 0)$p1
      gamma<-ifelse(is.na(gamma),yes = 10,no = gamma)
      g <- gamma * sqrt(n)
      d <- sqrt(sum((log(x) - mu)^2) + delta^2) / sqrt(n)
    }else{
      gamma <- (3 / sqrt(n))
      s2 <- ifelse(is.numeric(guess_s2), yes = guess_s2, no = s2)
      delta <- optimx::optimx(functional_delta, par = 1, lower = 0, method = "nlminb", n = n,
                              quant = quant, gamma = gamma, sigma2 = s2, lambda = 0)$p1
      delta<-ifelse(is.na(delta),yes = 1,no = delta)
      g <- gamma * sqrt(n)
      d <- sqrt(sum((log(x) - mu)^2) + delta^2) / sqrt(n)
    }
    par <- round(matrix(c(0, delta, gamma, NA, NA, l, d, g, mu, beta), ncol = 5, byrow = T), 3)
  }else{
    stop("The selected method must be 'weak_inf' or 'optimal'")
  }

  if(quant == 0.5){
    est <- SMNG_MGF(r = 1, mu = -mu, delta = d, gamma = g, lambda = l, beta = beta, inf_sum = T)/
      SMNG_MGF(r = 1, mu = - 2 * mu, delta = 2 * d, gamma = g / 2, lambda = l, beta = beta, inf_sum = T)
    message("The Bayes estimator under relative quadratic loss is employed\n")
  }else{
    est <- SMNG_MGF(r = 1, mu = mu, delta = d, gamma = g, lambda = l, beta = beta, inf_sum = T)
    message("The Bayes estimator under quadratic loss is employed\n")
  }
  var <- SMNG_MGF(r = 2, mu = mu, delta = d, gamma = g, lambda = l, beta = beta, inf_sum = T) -
    SMNG_MGF(r = 1, mu = mu, delta = d, gamma = g, lambda = l, beta = beta, inf_sum = T) ^ 2
  colnames(par) <- c("lambda", "delta", "gamma", "mu", "beta")
  rownames(par) <- c("prior", "posterior")

  if(CI == FALSE){
    xi_post=ghyp::ghyp.ad(mu = mu, delta = d, alpha = g,
                         lambda = l, beta = 0)
    log_par<-matrix(nrow = 2, ncol = 5)
    log_par[1,1]<-ghyp::ghyp.moment(xi_post, order = 1, central = F)
    log_par[1,2]<-ghyp::ghyp.moment(xi_post, order = 2, central = T)
    log_par[1,3:5]<-ghyp::qghyp(p = c(0.05, 0.5, 0.95), xi_post)
    log_par[2,1]<-ghyp::Egig(lambda = l, chi = d^2*n,psi = (g/sqrt(n))^2, func = "x")
    log_par[2,2]<-ghyp::Egig(lambda = l, chi = d^2*n,psi = (g/sqrt(n))^2, func = "var")
    log_par[2,3:5]<-ghyp::qgig(p = c(0.05, 0.5, 0.95), lambda = l, chi = d^2*n,psi = (g/sqrt(n))^2)
    colnames(log_par)<-c("Mean","Var","p=0.05","p=0.50","p=0.95")
    rownames(log_par)<-c("xi","sigma2")
    return(list(Quantile = quant, Parameters = par,
                LogN_Par_Post=log_par, Post_Estimates = as.matrix(data.frame(Mean=est, S.d.=sqrt(var)))))
  } else {
    if (!(method_CI == "exact" | method_CI == "simulation"))
      stop("method_CI must be between 'exact' or 'simulation'")
    if (type_CI == "two-sided") {
      if (method_CI == "exact") {
        low <- qlSMNG(lambda = l, gamma = g, delta = d, beta = beta, mu = mu, p = alpha_CI / 2, rel_tol = rel_tol_CI)
        up <- qlSMNG(lambda = l, gamma = g, delta = d, beta = beta, mu = mu, p = (1 - alpha_CI / 2), rel_tol = rel_tol_CI)
      } else {
        sample <- rlSMNG(n = nrep_CI, lambda = l, gamma = g, delta = d, beta = beta, mu = mu)
        low <- quantile(x = sample, probs = alpha_CI / 2)
        up <- quantile(x = sample, probs = (1 - alpha_CI / 2))
      }
    } else if (type_CI == "LCL") {
      if (method_CI == "exact") {
        low <- qlSMNG(lambda = l, gamma = g, delta = d, beta = beta, mu = mu, p = alpha_CI, rel_tol = rel_tol_CI)
      } else {
        sample <- rlSMNG(n = nrep_CI, lambda = l, gamma = g, delta = d, beta = beta, mu = mu)
        low <- quantile(x = sample, probs = alpha_CI)
      }
      up <- Inf
    } else if (type_CI == "UCL") {
      low <- 0
      if (method_CI == "exact") {
        up <- qlSMNG(lambda = l, gamma = g, delta = d, beta = beta, mu = mu, p = (1 - alpha_CI), rel_tol = rel_tol_CI)
      } else {
        sample <- rlSMNG(n = nrep_CI, lambda = l, gamma = g, delta = d, beta = beta, mu = mu)
        up <- quantile(x = sample, probs = (1 - alpha_CI))
      }
    } else {
      stop("type_CI must be 'two-sided', 'LCL' or 'UCL'")
    }
    xi_post=ghyp::ghyp.ad(mu = mu, delta = d, alpha = g,
                          lambda = l, beta = 0)
    log_par<-matrix(nrow = 2, ncol = 5)
    log_par[1,1]<-ghyp::ghyp.moment(xi_post, order = 1, central = F)
    log_par[1,2]<-ghyp::ghyp.moment(xi_post, order = 2, central = T)
    log_par[1,3:5]<-ghyp::qghyp(p = c(0.05, 0.5, 0.95), xi_post)
    log_par[2,1]<-ghyp::Egig(lambda = l, chi = d^2*n,psi = (g/sqrt(n))^2, func = "x")
    log_par[2,2]<-ghyp::Egig(lambda = l, chi = d^2*n,psi = (g/sqrt(n))^2, func = "var")
    log_par[2,3:5]<-ghyp::qgig(p = c(0.05, 0.5, 0.95), lambda = l, chi = d^2*n,psi = (g/sqrt(n))^2)
    colnames(log_par)<-c("Mean","Var","p=0.05","p=0.50","p=0.95")
    rownames(log_par)<-c("xi","sigma2")
    limits<- c(low, up)
    names(limits)<-c("Lower limit", "Upper limit")
    return(list(Quantile = quant, Parameters = par, LogN_Par_Post=log_par, Post_Estimates = as.matrix(data.frame(Mean=est, S.d.=sqrt(var))), Interval = limits))

  }

}







#'Bayesian estimate of the log-normal conditioned quantiles
#'
#'This function produces a point estimate for the log-normal distribution quantile of fixed level \code{quant}.
#'
#'@param y Vector of observations of the response variable.
#'@param X Design matrix.
#'@param Xtilde Covariate patterns of the units to estimate.
#'@param quant Number between 0 and 1 that indicates the quantile of interest.
#'@param method String that indicates the prior setting to adopt. Choosing \code{"weak_inf"}
#'a weakly informative prior setting is adopted, whereas selecting
#' \code{"optimal"} the hyperparameters are fixed trough a numerical optimization algorithm
#' aimed at minimizing the frequentist MSE.
#' @param guess_s2 Specification of a guess for the variance if available. If not, the sample estimate is used.
#' @param y_transf Logical. If \code{TRUE}, the \code{y} vector is assumed already log-transformed.
#' @param CI Logical. With the default choice \code{TRUE}, the posterior credibility interval is computed.
#'@param method_CI String that indicates if the limits should be computed through the logSMNG
#'quantile function \code{\link{qlSMNG}} (option \code{"exact"}, default), or by randomly generating
#'(\code{"simulation"}) using the function \code{\link{rlSMNG}}.
#' @param alpha_CI Level of credibility of the posterior interval.
#' @param type_CI String that indicates the type of interval to compute: \code{"two-sided"} (default),
#'\code{"UCL"} (i.e. Upper Credible Limit) for upper one-sided intervals  or \code{"LCL"} (i.e. Lower
#'Credible Limit) for lower one-sided intervals.
#'@param rel_tol_CI Level of relative tolerance required for the \code{integrate} procedure or for the infinite sum.
#' Default set to \code{1e-5}.
#'@param nrep Number of simulations for the C.I. in case of \code{method="simulation"} and for the posterior of the coefficients vector.
#'
#'@details
#' The function allows to carry out Bayesian inference for the conditional quantiles of a sample that is assumed log-normally distributed.
#' The design matrix containing the covariate patterns of the sampled units is \code{X}, whereas \code{Xtilde}
#' contains the covariate patterns of the unit to predict.
#'
#' The classical log-normal linear mixed model is assumed and the quantiles are estimated as:
#' \deqn{\theta_p(x)=exp(x^T\beta+\Phi^{-1}(p))}.
#'
#' A generalized inverse Gaussian prior is assumed for the variance in the log scale \eqn{\sigma^2}, whereas a
#' flat improper prior is assumed for the vector of coefficients \eqn{\beta}.
#'
#' Two alternative hyperparamters setting are implemented (choice controlled by the argument \code{method}): a weakly
#' informative proposal and an optimal one.
#'
#'
#'
#'@return
#' The function returns the prior parameters and their posterior values, summary statistics of the parameters \eqn{\beta} and \eqn{\sigma^2}, and the estimate of the specified quantile:
#' the posterior mean and variance are provided by default. Moreover the user can control the computation of posterior intervals.
#'
#'#'@source
#'
#'Gardini, A., C. Trivisano, and E. Fabrizi. \emph{Bayesian inference for quantiles of the log-normal distribution.} Biometrical Journal (2020).
#'
#'
#'
#' @export


LN_QuantReg <- function(y, X, Xtilde, quant, method = "weak_inf", guess_s2=NULL, y_transf=TRUE, CI = TRUE,
                        method_CI = "exact",alpha_CI = 0.05, type_CI = "two-sided",rel_tol_CI = 1e-5, nrep = 100000) {
  if(is.vector(y) == FALSE)
    stop("y must be a vector")
  if(is.matrix(Xtilde) == FALSE)
    stop("Xtilde must be a matrix")
  if(is.matrix(X) == FALSE)
    stop("X must be a matrix")
  if(length(y) != dim(X)[1])
    stop("y and X must have the same sample size")
  if(dim(Xtilde)[2] != dim(X)[2])
    stop("Xtilde and X must have the same number of variables")
  if (!(method == "weak_inf" | method == "optimal"))
    stop("method must be between 'weak_inf' or 'optimal'")
  n <- length(y)
  if (y_transf == FALSE){
    z <- log(y)
  } else {
    z <- y
  }
  l_mod <- lm(z ~ -1+X)
  mu <- as.numeric(Xtilde %*% l_mod$coefficients)
  RSS <- sum(l_mod$residuals ^ 2)
  nreg <- dim(X)[2]
  s2 <- RSS / (n - nreg)
  lev<-numeric(dim(Xtilde)[1])
  for(j in 1:dim(Xtilde)[1]){
    lev[j]<-Xtilde[j,]%*%solve(t(X)%*%X)%*%t(Xtilde)[,j]
  }
  h_M <- max(lev)
  h <- nreg / n
  beta <- qnorm(p = quant) / sqrt(h)

  if(method == "weak_inf"){
    l <- -(n / 2 - nreg / 2)
    gamma<-(3 * sqrt(h_M))
    g <- gamma / sqrt(h)
    delta <- 0.01
    d <- sqrt(RSS + delta ^ 2) * sqrt(h)
    par <- round(matrix(c(0, 0.01,gamma,
                          NA, l, d, g, beta), ncol = 4, byrow = T), 3)
  }else if(method=="optimal"){
    l<--(n/2-nreg / 2)
    if(quant==0.5){
      stop("The 'optimal' method cannot be used for the median: use the 'weak_inf' method")
    }else if(quant > 0.5){
      g0 <- (3 * sqrt(h_M))
      delta <- ifelse(is.numeric(guess_s2), yes = guess_s2, no = 1)
      s2 <- ifelse(is.numeric(guess_s2), yes = guess_s2, no = s2)
      gamma <- optimx::optimx(fn = functional_gamma, par = g0 + 1,lower = g0, method = "nlminb", n = n,
                              quant = quant, delta = delta, nreg = nreg, sigma2 = s2, lambda = 0)$p1
      if(is.na(gamma)){
        warning("Failed to find the optimal value for gamma because of the smallness of sigma2:
                substituted with 15")
        gamma <- 15
        }
      g <- gamma / sqrt(h)
      d <- sqrt(RSS + delta^2) * sqrt(h)
    }else{
      gamma <- (3 * sqrt(h_M))
      s2 <- ifelse(is.numeric(guess_s2), yes = guess_s2, no = s2)
      delta <- optimx::optimx(functional_delta, par = 1, lower = 0, method = "nlminb", n = n,
                              quant = quant, gamma = gamma, sigma2 = s2, nreg = nreg, lambda = 0)$p1
      if(is.na(delta)){
        warning("Failed to find the optimal value for delta because of the smallness of sigma2:
                substituted with 15")
        delta <- 15

      }
      g <- gamma / sqrt(h)
      d <- sqrt(RSS + delta^2) * sqrt(h)
    }
    par <- round(matrix(c(0, delta, gamma, NA, l, d, g, beta), ncol = 4, byrow = T), 3)
  }else{
    stop("The selected method must be 'weak_inf' or 'optimal'")
  }
  n_pred<-dim(Xtilde)[1]
  est<-numeric(n_pred)
  var<-numeric(n_pred)

  if(quant == 0.5){
    for(i in 1:n_pred){
    est[i] <- SMNG_MGF(r = 1, mu = -mu[i], delta = d, gamma = g, lambda = l, beta = beta, inf_sum = T)/
      SMNG_MGF(r = 1, mu = - 2 * mu[i], delta = 2 * d, gamma = g / 2, lambda = l, beta = beta, inf_sum = T)
    }
    message("The Bayes estimator under relative quadratic loss is employed\n")
  }else{
    for(i in 1:n_pred){
    est[i] <- SMNG_MGF(r = 1, mu = mu[i], delta = d, gamma = g, lambda = l, beta = beta, inf_sum = T)
    }
    message("The Bayes estimator under quadratic loss is employed\n")
  }
  for(i in 1:n_pred){
  var[i] <- SMNG_MGF(r = 2, mu = mu[i], delta = d, gamma = g, lambda = l, beta = beta, inf_sum = T) -
    SMNG_MGF(r = 1, mu = mu[i], delta = d, gamma = g, lambda = l, beta = beta, inf_sum = T) ^ 2
  }
  colnames(par) <- c("lambda", "delta", "gamma", "beta")
  rownames(par) <- c("prior", "posterior")



  beta_post <- matrix(nrow = nreg, ncol = 7)
  colnames(beta_post) <- c("Mean", "S.d.", "q2.5", "q25", "q50", "q75", "q97.5")
  rownames(beta_post)<-colnames(X)
  V<-solve(t(X)%*%X)

  for(j in 1:nreg){
    gen <- suppressWarnings(ghyp::rghyp(n = nrep, object = ghyp::ghyp.ad(lambda = l,
                                                        alpha = gamma / sqrt(diag(V)[j]), delta = sqrt(RSS + delta^2)*sqrt(diag(V)[j]),
                                                        mu =  as.numeric(l_mod$coefficients[j]))))
    beta_post[j,1]<-mean(gen)
    beta_post[j,2]<-sd(gen)
    beta_post[j,3:7]<- quantile(gen, probs = c(0.025, 0.25, 0.5, 0.75, 0.975))
    gen<-NULL
}

  if(CI == FALSE){
    sigma2<-matrix(nrow = 1, ncol = 5)
    sigma2[1,1]<-ghyp::Egig(lambda = l, chi = d^2/h,psi = (g*sqrt(h))^2, func = "x")
    sigma2[1,2]<-ghyp::Egig(lambda = l, chi = d^2/h,psi = (g*sqrt(h))^2, func = "var")
    sigma2[1,3:5]<-ghyp::qgig(p = c(0.05, 0.5, 0.95), lambda = l, chi = d^2/h,psi = (g*sqrt(h))^2)
    colnames(sigma2)<-c("Mean","Var","q5","q50","q95")
    rownames(sigma2)<-c("sigma2")
    return(list(Quantile = quant, Parameters = par, Leverage = h, Sigma2=sigma2, Coefficients = beta_post, Post_Estimates = as.matrix(data.frame(Mean=est, S.d.=sqrt(var)))))
  } else {
    low<-numeric(n_pred)
    up<-numeric(n_pred)
    if (!(method_CI == "exact" | method_CI == "simulation"))
      stop("method_CI must be between 'exact' or 'simulation'")
    if (type_CI == "two-sided") {
      if (method_CI == "exact") {
        for(i in 1:n_pred){
        low[i] <- qlSMNG(lambda = l, gamma = g, delta = d, beta = beta, mu = mu[i], p = alpha_CI / 2, rel_tol = rel_tol_CI)
        up[i] <- qlSMNG(lambda = l, gamma = g, delta = d, beta = beta, mu = mu[i], p = (1 - alpha_CI / 2), rel_tol = rel_tol_CI)
        }
        } else {
          for(i in 1:n_pred){
        sample <- rlSMNG(n = nrep, lambda = l, gamma = g, delta = d, beta = beta, mu = mu[i])
        low[i] <- quantile(x = sample, probs = alpha_CI / 2)
        up[i] <- quantile(x = sample, probs = (1 - alpha_CI / 2))
          }
      }
    } else if (type_CI == "LCL") {
      for(i in 1:n_pred){
      if (method_CI == "exact") {
        low[i] <- qlSMNG(lambda = l, gamma = g, delta = d, beta = beta, mu = mu[i], p = alpha_CI, rel_tol = rel_tol_CI)
      } else {
        sample <- rlSMNG(n = nrep, lambda = l, gamma = g, delta = d, beta = beta, mu = mu[i])
        low[i] <- quantile(x = sample, probs = alpha_CI)
      }
      up[i] <- Inf
      }
    } else if (type_CI == "UCL") {
      for(i in 1:n_pred){
      low[i] <- 0
      if (method_CI == "exact") {
        up[i] <- qlSMNG(lambda = l, gamma = g, delta = d, beta = beta, mu = mu[i], p = (1 - alpha_CI), rel_tol = rel_tol_CI)
      } else {
        sample <- rlSMNG(n = nrep, lambda = l, gamma = g, delta = d, beta = beta, mu = mu[i])
        up[i] <- quantile(x = sample, probs = (1 - alpha_CI))
      }
      }
    } else {
      stop("type_CI must be 'two-sided', 'LCL' or 'UCL'")
    }
    limits<- cbind(low, up)
    names(limits)<-c("Lower limit", "Upper limit")
    sigma2<-matrix(nrow = 1, ncol = 5)
    sigma2[1,1]<-ghyp::Egig(lambda = l, chi = d^2/h,psi = (g*sqrt(h))^2, func = "x")
    sigma2[1,2]<-ghyp::Egig(lambda = l, chi = d^2/h,psi = (g*sqrt(h))^2, func = "var")
    sigma2[1,3:5]<-ghyp::qgig(p = c(0.05, 0.5, 0.95), lambda = l, chi = d^2/h,psi = (g*sqrt(h))^2)
    colnames(sigma2)<-c("Mean","Var","q5","q50","q95")
    rownames(sigma2)<-c("sigma2")
    return(list(Quantile = quant, Parameters = par, Sigma2=sigma2, Coefficients = beta_post,
                Post_Estimates = as.matrix(data.frame(Mean=est, S.d.=sqrt(var))), Interval = limits[,]))
  }

}

