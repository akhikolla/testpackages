

#---------------------------------------------
#' ZiBell regression model
#' @aliases zibellreg
#' @export
#' @description Fits the Bell regression model to overdispersed count data.
#' @param formula an object of class "formula" (or one that can be coerced to that class): a symbolic description of the model to be fitted.
#' @param data an optional data frame, list or environment (or object coercible by as.data.frame to a data frame) containing the variables in the model. If not found in data, the variables are taken from environment(formula), typically the environment from which ypbp is called.
#' @param approach approach to be used to fit the model (mle: maximum likelihood; bayes: Bayesian approach).
#' @param hessian hessian logical; If TRUE (default), the hessian matrix is returned when approach="mle".
#' @param hyperpars a list containing the hyperparameters associated with the prior distribution of the regression coefficients; if not specified then default choice is hyperpars = c(mu_psi = 0, sigma_psi = 10, mu_beta = 0, sigma_beta = 10).
#' @param ... further arguments passed to either `rstan::optimizing` or `rstan::sampling`.
#' @return zibellreg returns an object of class "zibellreg" containing the fitted model.
#'
#' @examples
#' \donttest{
#' # ML approach:
#' mle <- zibellreg(cells ~ smoker+gender|smoker+gender, data = cells, approach = "mle")
#' summary(mle)
#'
#' # Bayesian approach:
#' bayes <- zibellreg(cells ~ 1|smoker+gender, data = cells, approach = "bayes")
#' summary(bayes)
#' }
#'
zibellreg<- function(formula, data, approach = c("mle", "bayes"), hessian = TRUE,
                   hyperpars = list(mu_psi=0, sigma_psi=10, mu_beta=0, sigma_beta=10), ...){
  approach <- match.arg(approach)
  formula <- Formula::Formula(formula)
  mf <- stats::model.frame(formula=formula, data=data)
  Terms <- stats::terms(mf)
  Z <- stats::model.matrix(formula, data = mf, rhs = 1)
  X <- stats::model.matrix(formula, data = mf, rhs = 2)
  Xlabels <- colnames(X)
  Zlabels <- colnames(Z)
  y <- stats::model.response(mf)
  n <- nrow(X)
  p <- ncol(X)
  q <- ncol(Z)

  if(p>1){
    if(match("(Intercept)", Xlabels)==1){
      X_std <- scale(X[,-1])
      x_mean <- array(c(0, attr(X_std, "scaled:center")))
      x_sd <- array(c(1, attr(X_std, "scaled:scale")))
      X_std <- cbind(1, X_std)
      Delta_x <- diag(1/x_sd)
      Delta_x[1,] <- Delta_x[1,] -  x_mean/x_sd
    }else{
      X_std <- scale(X)
      x_mean <- array(attr(X_std, "scaled:center"))
      x_sd <- array(attr(X_std, "scaled:scale"))
      Delta_x <- diag(1/x_sd)
    }
  }else{
    X_std <- X
    x_mean <- array(0)
    x_sd <- array(1)
  }


  if(q>1){
    if(match("(Intercept)", Zlabels)==1){
      Z_std <- scale(Z[,-1])
      z_mean <- array(c(0, attr(Z_std, "scaled:center")))
      z_sd <- array(c(1, attr(Z_std, "scaled:scale")))
      Z_std <- cbind(1, Z_std)
      Delta_z <- diag(1/z_sd)
      Delta_z[1,] <- Delta_z[1,] -  z_mean/z_sd
    }else{
      Z_std <- scale(Z)
      z_mean <- array(attr(Z_std, "scaled:center"))
      z_sd <- array(attr(Z_std, "scaled:scale"))
      Delta_z <- diag(1/z_sd)
    }
  }else{
      Z_std <- Z
      z_mean <- array(0)
      z_sd <- array(1)
    }


  stan_data <- list(y=y, X=X_std, Z=Z_std, n=n, p=p, q=q, x_mean=x_mean, x_sd=x_sd, z_mean=z_mean, z_sd=z_sd,
                    mu_beta = hyperpars$mu_beta, sigma_beta=hyperpars$sigma_beta,
                    mu_psi = hyperpars$mu_psi, sigma_psi=hyperpars$sigma_psi,
                    approach=0)

  if(approach=="mle"){
    fit <- rstan::optimizing(stanmodels$zibellreg, hessian=hessian,
                             data=stan_data, verbose=FALSE, ...)
    if(hessian==TRUE){
      fit$hessian <- - fit$hessian
    }
    fit$par <- fit$theta_tilde[-(1:(p+q))]
    AIC <- -2*fit$value + 2*(p+q)
    fit <- list(fit=fit, logLik = fit$value, AIC = AIC, Delta = magic::adiag(Delta_z, Delta_x))
  }else{
    stan_data$approach <- 1
    fit <- rstan::sampling(stanmodels$zibellreg, data=stan_data, verbose=FALSE, ...)
    fit <- list(fit=fit)
  }


  fit$n <- n
  fit$p <- p
  fit$q <- q
  # fit$x_mean <- x_mean
  # fit$x_sd <- x_sd
  # fit$z_mean <- z_mean
  # fit$z_sd <- z_sd
  # fit$v_sd <- c(z_sd, x_sd)


  fit$call <- match.call()
  fit$formula <- stats::formula(Terms)
  fit$terms <- stats::terms.formula(formula)
  fit$labels1 <- Zlabels
  fit$labels2 <- Xlabels
  fit$approach <- approach
  class(fit) <- "zibellreg"
  return(fit)
}


