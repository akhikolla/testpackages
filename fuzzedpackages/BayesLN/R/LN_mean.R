#' Bayesian Estimate of the Log-normal Mean
#'
#' This function produces a Bayesian estimate of the log-normal mean, assuming a GIG prior for the variance and an
#' improper flat prior for the mean in the log scale.
#'
#' @param x Vector containing the sample.
#' @param method  String that indicates the prior setting to adopt. Choosing \code{"weak_inf"} a weakly informative prior setting is adopted, whereas selecting
#' \code{"optimal"} the hyperparameters are aimed at minimizing the frequentist MSE.
#' @param x_transf Logical. If \code{TRUE}, the \code{x} vector is assumed already log-transformed.
#' @param CI Logical. With the default choice \code{TRUE}, the posterior credibility interval is computed.
#' @param alpha_CI Level of alpha that determines the credibility (1-\code{alpha_CI}) of the posterior interval.
#' @param type_CI String that indicates the type of interval to compute: \code{"two-sided"} (default),
#'\code{"UCL"} (i.e. Upper Credible Limit) for upper one-sided intervals or \code{"LCL"} (i.e. Lower
#'Credible Limit) for lower one-sided intervals.
#'@param nrep Number of simulations for the computation of the credible intervals.
#'
#'@details Summarizing the posterior mean of the log-normal expectation might be delicate since several
#'common priors used for the variance do not produces posteriors with finite moments. The proposal by Fabrizi and Trivisano (2012) of adopting a generalized inverse Gaussian (GIG)
#'prior for the variance in the log scale \eqn{\sigma^2} has been implemented. Moreover, they discussed how to specify the hyperparameters according to two different aims.
#'
#'Firstly, a weakly informative
#'prior allowed to produce posterior credible intervals with good frequentist properties, whereas a prior aimed at minimizing the point estimator
#'MSE was proposed too. The choice between the two priors can be made through the argument \code{method}.
#'
#'The point estimates are exact values, whereas the credible intervals are provided through simulations from the posterior distribution.
#'
#'
#'@return The function returns a list which includes the prior and posterior parameters, the point estimate of the log-normal mean that consists in the mean of the posterior
#'distribution of the functional \eqn{\exp\{\mu+\sigma^2/2\}} and the posterior variance.
#'
#'@source Fabrizi, E., & Trivisano, C. \emph{Bayesian estimation of log-normal means with finite quadratic expected loss}. Bayesian Analysis, 7(4), 975-996. (2012).
#'
#'@examples
#' # Load data
#' data("NCBC")
#' # Optimal point estimator
#' LN_Mean(x = NCBC$al, x_transf = FALSE, method = "optimal", CI = FALSE)
#' # Weakly informative prior and interval estimation
#' LN_Mean(x = NCBC$al, x_transf = FALSE, type_CI = "UCL")
#'
#' @export


LN_Mean <- function (x, method = "weak_inf", x_transf = TRUE, CI = TRUE,
                     alpha_CI = 0.05, type_CI = "two-sided", nrep=1e5) {
  if (!(method != "optimal" | method != "weak_inf"))
    stop ("method can be 'optimal' or ' weak_inf'")
  if (x_transf == FALSE){
    x <- log(x)
  } else {
    x <- x
  }
  mu <- mean(x)
  s2 <- var(x)
  n <- length(x)
  g <- sqrt(9/n+3)
  d <- 0.01
  l_opt <- (n - 3) / 2 - (n ^ 2 - 1) / (2 * n - 6)
  l <- ifelse(test = method == "optimal", yes =  l_opt, no = 0)
  beta <- n * 0.5
  alpha <- sqrt(n * (g ^ 2 + n * 0.5 ^ 2))
  delta_bar <- sqrt((1 / n) * (s2 * (n - 1) + d ^ 2))
  l_bar <- l - n / 2 + 0.5

  est <- GH_MGF(r=1, mu = mu, delta = delta_bar, alpha = alpha,
                lambda = l_bar, beta = beta)

  var<- GH_MGF(r=2, mu = mu, delta = delta_bar, alpha = alpha,
              lambda = l_bar, beta = beta) - est^2
  prior.par <- round(c(l, d, g), 3)
  names(prior.par) <- c("lambda", "delta", "gamma")

  posterior.par <- round(c(l_bar, alpha, delta_bar, beta, mu),3)
  names(posterior.par) <- c("lambda", "alpha", "delta", "beta", "mu")

  xi_post=ghyp::ghyp.ad(mu = mu, delta = delta_bar, alpha = g * sqrt(n),
                        lambda = l_bar, beta = 0)
  log_par<-matrix(nrow = 2, ncol = 5)
  log_par[1,1]<-ghyp::ghyp.moment(xi_post, order = 1, central = F)
  log_par[1,2]<-ghyp::ghyp.moment(xi_post, order = 2, central = T)
  log_par[1,3:5]<-ghyp::qghyp(p = c(0.05, 0.5, 0.95), xi_post)
  log_par[2,1]<-ghyp::Egig(lambda = l_bar, chi = delta_bar^2*n,psi = g^2, func = "x")
  log_par[2,2]<-ghyp::Egig(lambda = l_bar, chi = delta_bar^2*n,psi = g^2, func = "var")
  log_par[2,3:5]<-ghyp::qgig(p = c(0.05, 0.5, 0.95), lambda = l_bar, chi = delta_bar^2*n,psi = g^2)
  colnames(log_par)<-c("Mean","Var","p=0.05","p=0.50","p=0.95")
  rownames(log_par)<-c("xi","sigma2")



  if(CI == FALSE){
      return(list(Prior_Parameters = prior.par,Posterior_Parameters = posterior.par, LogN_Par_Post=log_par,
                  Post_Estimates = as.matrix(data.frame(Mean=est, S.d.=sqrt(var)))))
  } else {
    if (type_CI == "two-sided") {
      sample <- suppressWarnings(exp(ghyp::rghyp(n = nrep,object=ghyp::ghyp.ad(mu = mu, delta = delta_bar, alpha = alpha,
                                                              lambda = l_bar, beta = beta))))
      low <- quantile(x = sample, probs = alpha_CI / 2)
      up <- quantile(x = sample, probs = (1 - alpha_CI / 2))
    } else if (type_CI == "LCL") {
      sample <- suppressWarnings(exp(ghyp::rghyp(n = nrep, object = ghyp::ghyp.ad(mu = mu, delta = delta_bar, alpha = alpha,
                                                              lambda = l_bar, beta = beta))))
      low <- quantile(x = sample, probs = alpha_CI)
      up <- Inf
    } else if (type_CI == "UCL") {
      low <- 0
      sample <- suppressWarnings(exp(ghyp::rghyp(n = nrep, object=ghyp::ghyp.ad(mu = mu, delta = delta_bar, alpha = alpha,
                                                              lambda = l_bar, beta = beta))))
      up <- quantile(x = sample, probs = (1 - alpha_CI))
    } else {
      stop("type must be 'two-sided', 'LCL' or 'UCL'")
    }
     limits<- c(low, up)
     names(limits)<-c("Lower limit", "Upper limit")
    return(list(Prior_Parameters = prior.par, Posterior_Parameters = posterior.par, LogN_Par_Post = log_par,
                Post_Estimates = as.matrix(data.frame(Mean=est, S.d.=sqrt(var))), Interval = limits))

  }
}





#' Bayesian Estimate of the conditional Log-normal Mean
#'
#' This function produces a bayesian estimate of the conditional log-normal mean assuming a GIG prior for the variance and an
#' improper prior for the regression coefficients of the linear regression in the log scale.
#'
#' @param y Vector of observations of the response variable.
#' @param X Design matrix.
#' @param Xtilde Matrix of covariate patterns for which an estimate is required.
#' @param method  String that indicates the prior setting to adopt. Choosing \code{"weak_inf"} a weakly
#' informative prior setting is adopted, whereas selecting
#' \code{"optimal"} the hyperparameters are aimed at minimizing the frequentist MSE.
#' @param y_transf Logical. If \code{TRUE}, the \code{y} vector is already assumed as log-transformed.
#' @param h Leverage. With the default option \code{NULL}, the average leverage is used.
#' @param CI Logical. With the default choice \code{TRUE}, the posterior credibility interval is computed.
#' @param alpha_CI Level of alpha that determines the credibility (1-\code{alpha_CI}) of the posterior interval.
#' @param type_CI String that indicates the type of interval to compute: \code{"two-sided"} (default),
#'\code{"UCL"} (i.e. Upper Credible Limit) for upper one-sided intervals  or \code{"LCL"} (i.e. Lower
#'Credible Limit) for lower one-sided intervals.
#' @param nrep Number of simulations.
#'
#'@details In this function the same procedure as \link{LN_Mean} is implemented allowing for the inclusion of covariates.
#'Bayesian point and interval estimates for the response variabile in the original scale are provided considering the model:
#'\eqn{log(y_i)=X\beta}.
#'
#'@return The function returns a list including the prior and posterior parameters, the point estimate of the log-normal mean conditioned with respect to the covariate
#'points included in \code{Xtilde}. It consists of the mean of the posterior
#'distribution for the functional \eqn{\exp\{\tilde{x}_i^T\beta+\sigma^2/2\}} and the posterior variance.
#'
#'@source Fabrizi, E., & Trivisano, C. \emph{Bayesian Conditional Mean Estimation in Log-Normal Linear Regression Models with Finite
#'Quadratic Expected Loss.} Scandinavian Journal of Statistics, 43(4), 1064-1077. (2016).
#'
#'@examples
#' library(BayesLN)
#' data("fatigue")
#'
#' # Design matrices
#' Xtot <- cbind(1, log(fatigue$stress), log(fatigue$stress)^2)
#' X <- Xtot[-c(1,13,22),]
#' y <- fatigue$cycle[-c(1,13,22)]
#' Xtilde <- Xtot[c(1,13,22),]
#' #Estimation
#' LN_MeanReg(y = y,
#'            X = X, Xtilde = Xtilde,
#'            method = "weak_inf", y_transf = FALSE)
#'
#'
#'@export


LN_MeanReg <- function(y, X, Xtilde, method = "weak_inf", y_transf=TRUE, h=NULL, CI = TRUE,
                       alpha_CI = 0.05, type_CI = "two-sided", nrep= 1e5) {
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

  h <- nreg / n
  m<-max(diag(Xtilde%*%solve(t(X)%*%X)%*%t(Xtilde)))
  g <- sqrt(3+9*m)
  d <- 0.01
  l_opt <- (n - nreg -2) / 2 - ((h + 1) * (n - nreg)) / (2 * (1 - 3 * h))
  l <- ifelse(test = method == "optimal", yes =  l_opt, no = 0)
  delta_bar <- sqrt(h * (RSS + d^2))
  l_bar <- l - (n - nreg) / 2
  beta <- 1 / (2 * h)
  alpha <- (1 / h) * sqrt(g ^ 2 + 1 / (4 * h))

  n_pred<-dim(Xtilde)[1]
  est<-numeric(n_pred)
  var<-numeric(n_pred)

  for(i in 1:n_pred){
    est[i] <- GH_MGF(r=1, mu = mu[i], delta = delta_bar, alpha = alpha,
                lambda = l_bar, beta = beta)

  var[i]<- GH_MGF(r=2, mu = mu[i], delta = delta_bar, alpha = alpha,
               lambda = l_bar, beta = beta) - est[i]^2
}
  prior.par <- round(c(l, d, g), 3)
  names(prior.par) <- c("lambda", "delta", "gamma")

  posterior.par <- round(cbind(l_bar, alpha, delta_bar, beta, mu),3)
  colnames(posterior.par) <- c("lambda", "alpha", "delta", "beta", "mu")

  beta_post <- matrix(nrow = nreg, ncol = 7)
  colnames(beta_post) <- c("Mean", "S.d.", "q2.5", "q25", "q50", "q75", "q97.5")
  rownames(beta_post)<-colnames(X)
  V<-solve(t(X)%*%X)

  for(j in 1:nreg){
    gen <- suppressWarnings(ghyp::rghyp(n = nrep, object = ghyp::ghyp.ad(lambda = l_bar,
                                                                         alpha = g / sqrt(diag(V)[j]), delta = sqrt(RSS + d^2)*sqrt(diag(V)[j]),
                                                                         mu =  as.numeric(l_mod$coefficients[j]))))
    beta_post[j,1]<-mean(gen)
    beta_post[j,2]<-sd(gen)
    beta_post[j,3:7]<- quantile(gen, probs = c(0.025, 0.25, 0.5, 0.75, 0.975))
    gen<-NULL
  }
  sigma2<-matrix(nrow = 1, ncol = 5)
  sigma2[1,1]<-ghyp::Egig(lambda = l_bar, chi = delta_bar^2/h,psi = g^2, func = "x")
  sigma2[1,2]<-ghyp::Egig(lambda = l_bar, chi = delta_bar^2/h, psi = g^2, func = "var")
  sigma2[1,3:5]<-ghyp::qgig(p = c(0.05, 0.5, 0.95), lambda = l_bar, chi = delta_bar^2/h, psi = g^2)
  colnames(sigma2)<-c("Mean","Var","q5","q50","q95")
  rownames(sigma2)<-c("sigma2")

  if(CI == FALSE){
    return(list(Prior_Parameters = prior.par, Posterior_Parameters = posterior.par, Leverage = h, Sigma2 = sigma2,
                Coefficients = beta_post, Post_Estimates = as.matrix(data.frame(Mean=est, S.d.=sqrt(var)))))
  } else {
    low<-numeric(n_pred)
    up<-numeric(n_pred)
for(i in 1:n_pred){
    if (type_CI == "two-sided") {
      sample <-  suppressWarnings(exp(ghyp::rghyp(n = nrep,object=ghyp::ghyp.ad(mu = mu[i], delta = delta_bar, alpha = alpha,
                                                              lambda = l_bar, beta = beta))))
      low[i] <- quantile(x = sample, probs = alpha_CI / 2)
      up[i] <- quantile(x = sample, probs = (1 - alpha_CI / 2))
    } else if (type_CI == "LCL") {
      sample <-  suppressWarnings(exp(ghyp::rghyp(n = nrep, object = ghyp::ghyp.ad(mu = mu[i], delta = delta_bar, alpha = alpha,
                                                                 lambda = l_bar, beta = beta))))
      low[i] <- quantile(x = sample, probs = alpha_CI)
      up[i] <- Inf
    } else if (type_CI == "UCL") {
      low[i] <- 0
      sample <- suppressWarnings(exp(ghyp::rghyp(n = nrep, object=ghyp::ghyp.ad(mu = mu[i], delta = delta_bar, alpha = alpha,
                                                               lambda = l_bar, beta = beta))))
      up[i] <- quantile(x = sample, probs = (1 - alpha_CI))
    } else {
      stop("type must be 'two-sided', 'LCL' or 'UCL'")
    }
  }
    limits<- cbind(low, up)
    names(limits)<-c("Lower limit", "Upper limit")
    return(list(Prior_Parameters = prior.par, Posterior_Parameters = posterior.par, Sigma2 = sigma2, Coefficients = beta_post,
                Post_Estimates = as.matrix(data.frame(Mean=est, S.d.=sqrt(var))) ,Interval = limits[,]))

  }
}


