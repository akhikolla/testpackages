#' Binary Gaussian Process (with/without time-series)
#'
#' @description The function fits Gaussian process models with binary response. The models can also include time-series term if a time-series binary response is observed. The estimation methods are based on PQL/PQPL and REML (for mean function and correlation parameter, respectively).
#'
#' @param X a design matrix with dimension \code{n} by \code{d}.
#' @param Y a response matrix with dimension \code{n} by \code{T}. The values in the matrix are 0 or 1. If no time-series, \code{T = 1}. If time-series is included, i.e., \code{T > 1}, the first column is the response vector of time 1, and second column is the response vector of time 2, and so on.
#' @param R a positive integer specifying the order of autoregression only if time-series is included. The default is 1.
#' @param L a positive integer specifying the order of interaction between \code{X} and previous \code{Y} only if time-series is included. The value couldn't nbe larger than R. The default is 1.
#' @param corr a list of parameters for the specifing the correlation to be used. Either exponential correlation function or Matern correlation function can be used. See \code{\link[GPfit]{corr_matrix}} for details.
#' @param nugget a positive value to use for the nugget. If \code{NULL}, a nugget will be estimated. The default is 1e-10.
#' @param constantMean logical. \code{TRUE} indicates that the Gaussian process will have a constant mean, plus time-series structure if \code{R} or \code{T} is greater than one; otherwise the mean function will be a linear regression in X, plus an intercept term and time-series structure.
#' @param orthogonalGP logical. \code{TRUE} indicates that the orthogonal Gaussian process will be used. Only available when \code{corr} is \code{list(type = "exponential", power = 2)}.
#' @param converge.tol convergence tolerance. It converges when relative difference with respect to \eqn{\eta_t} is smaller than the tolerance. The default is 1e-10.
#' @param maxit a positive integer specifying the maximum number of iterations for estimation to be performed before the estimates are convergent. The default is 15.
#' @param maxit.PQPL a positive integer specifying the maximum number of iterations for PQL/PQPL estimation to be performed before the estimates are convergent. The default is 20.
#' @param maxit.REML a positive integer specifying the maximum number of iterations in \code{lbfgs} for REML estimation to be performed before the estimates are convergent. The default is 100.
#' @param xtol_rel a postive value specifying the convergence tolerance for \code{lbfgs}. The default is 1e-10.
#' @param standardize logical. If \code{TRUE}, each column of X will be standardized into \code{[0,1]}. The default is \code{FALSE}.
#' @param verbose logical. If \code{TRUE}, additional diagnostics are printed. The default is \code{TRUE}.
#'
#' @details Consider the model \deqn{logit(p_t(x))=\eta_t(x)=\sum^R_{r=1}\varphi_ry_{t-r}(\mathbf{x})+\alpha_0+\mathbf{x}'\boldsymbol{\alpha}+\sum^L_{l=1}\boldsymbol{\gamma}_l\textbf{x}y_{t-l}(\mathbf{x})+Z_t(\mathbf{x}),} where \eqn{p_t(x)=Pr(y_t(x)=1)} and \eqn{Z_t(\cdot)} is a Gaussian process with zero mean, variance \eqn{\sigma^2}, and correlation function \eqn{R_{\boldsymbol{\theta}}(\cdot,\cdot)} with unknown correlation parameters \eqn{\boldsymbol{\theta}}. The power exponential correlation function is used and defined by \deqn{R_{\boldsymbol{\theta}}(\mathbf{x}_i,\mathbf{x}_j)=\exp\{-\sum^{d}_{l=1}\frac{(x_{il}-x_{jl})^p}{\theta_l} \},} where \eqn{p} is given by \code{power}. The parameters in the mean functions including \eqn{\varphi_r,\alpha_0,\boldsymbol{\alpha},\boldsymbol{\gamma}_l} are estimated by PQL/PQPL, depending on whether the mean function has the time-series structure. The parameters in the Gaussian process including \eqn{\boldsymbol{\theta}} and \eqn{\sigma^2} are estimated by REML. If \code{orthogonalGP} is \code{TRUE}, then orthogonal covariance function (Plumlee and Joseph, 2016) is employed. More details can be seen in Sung et al. (unpublished).
#'
#' @return
#' \item{alpha_hat}{a vector of coefficient estimates for the linear relationship with X.}
#' \item{varphi_hat}{a vector of autoregression coefficient estimates.}
#' \item{gamma_hat}{a vector of the interaction effect estimates.}
#' \item{theta_hat}{a vector of correlation parameters.}
#' \item{sigma_hat}{a value of sigma estimate (standard deviation).}
#' \item{nugget_hat}{if \code{nugget} is \code{NULL}, then return a nugget estimate, otherwise return \code{nugget}.}
#' \item{orthogonalGP}{\code{orthogonalGP}.}
#' \item{X}{data \code{X}.}
#' \item{Y}{data \code{Y}.}
#' \item{R}{order of autoregression.}
#' \item{L}{order of interaction between X and previous Y.}
#' \item{corr}{a list of parameters for the specifing the correlation to be used.}
#' \item{Model.mat}{model matrix.}
#' \item{eta_hat}{eta_hat.}
#' \item{standardize}{\code{standardize}.}
#' \item{mean.x}{mean of each column of \code{X}. Only available when \code{standardize=TRUE}, otherwise \code{mean.x=NULL}.}
#' \item{scale.x}{\code{max(x)-min(x)} of each column of \code{X}. Only available when \code{standardize=TRUE}, otherwise \code{scale.x=NULL}.}
#'
#' @seealso \code{\link{predict.binaryGP}} for prediction of the binary Gaussian process, \code{\link{print.binaryGP}} for the list of estimation results, and \code{\link{summary.binaryGP}} for summary of significance results.
#' @author Chih-Li Sung <iamdfchile@gmail.com>
#'
#' @import nloptr
#' @import lhs
#' @import GPfit
#' @import utils
#' @import stats
#' @import graphics
#' @import methods
#' @importFrom Rcpp evalCpp
#'
#'
#' @examples
#' library(binaryGP)
#'
#' #####      Testing function: cos(x1 + x2) * exp(x1*x2) with TT sequences      #####
#' #####   Thanks to Sonja Surjanovic and Derek Bingham, Simon Fraser University #####
#' test_function <- function(X, TT)
#' {
#'   x1 <- X[,1]
#'   x2 <- X[,2]
#'
#'   eta_1 <- cos(x1 + x2) * exp(x1*x2)
#'
#'   p_1 <- exp(eta_1)/(1+exp(eta_1))
#'   y_1 <- rep(NA, length(p_1))
#'   for(i in 1:length(p_1)) y_1[i] <- rbinom(1,1,p_1[i])
#'   Y <- y_1
#'   P <- p_1
#'   if(TT > 1){
#'     for(tt in 2:TT){
#'       eta_2 <- 0.3 * y_1 + eta_1
#'       p_2 <- exp(eta_2)/(1+exp(eta_2))
#'       y_2 <- rep(NA, length(p_2))
#'       for(i in 1:length(p_2)) y_2[i] <- rbinom(1,1,p_2[i])
#'       Y <- cbind(Y, y_2)
#'       P <- cbind(P, p_2)
#'       y_1 <- y_2
#'     }
#'   }
#'
#'   return(list(Y = Y, P = P))
#' }
#'
#' set.seed(1)
#' n <- 30
#' n.test <- 10
#' d <- 2
#' X <- matrix(runif(d * n), ncol = d)
#'
#' ##### without time-series #####
#' Y <- test_function(X, 1)$Y  ## Y is a vector
#'
#' binaryGP.model <- binaryGP_fit(X = X, Y = Y)
#' print(binaryGP.model)   # print estimation results
#' summary(binaryGP.model) # significance results
#' \donttest{
#' ##### with time-series, lag 1 #####
#' Y <- test_function(X, 10)$Y  ## Y is a matrix with 10 columns
#'
#' binaryGP.model <- binaryGP_fit(X = X, Y = Y, R = 1)
#' print(binaryGP.model)   # print estimation results
#' summary(binaryGP.model) # significance results
#' }
#'
#' @useDynLib binaryGP, .registration = TRUE
#' @export

#@exportPattern "^[[:alpha:]]+"

binaryGP_fit <- function(X, Y, R = 0, L = 0, corr = list(type = "exponential", power = 2), nugget = 1e-10,
                         constantMean = FALSE, orthogonalGP = FALSE, converge.tol = 1e-10, maxit = 15, maxit.PQPL = 20, maxit.REML = 100, xtol_rel = 1e-10,
                         standardize = FALSE, verbose = TRUE){

  ###################       Setting       ###################
  #if (L > R) stop("L should be smaller than or equal to R")
  if(corr$type != "exponential" & corr$type != "matern") stop("corr$type can only be either exponential or matern.")
  if(corr$type == "exponential") if(corr$power > 2 | corr$power < 1) stop("corr$power can only be between 1 and 2.")
  T_val <- ncol(Y)
  if(is.null(T_val)) T_val <- 1       ### without time-series
  if(T_val == 1 & (R > 0 | L > 0) ){
    warning("R = 0 and L = 0 when T = 1.")
    R <- 0
    L <- 0
  }
  if(T_val > 1){
    if (R >= T_val) stop("R >= T")
    if (L >= T_val) stop("L >= T")
  }

  if (length(unique(as.numeric(Y))) != 2) stop("Number of response types is greater/less than 2.")
  if(any(Y != 1 & Y != 0)) stop("The value of matrix Y should be 0 or 1.")
  if(orthogonalGP & !(corr$type == "exponential" & corr$power == 2))  {
    warning("orthogonalGP = TRUE only can be used for power exponential correlation function with power = 2, so orthogonalGP = FALSE \n")
    orthogonalGP <- FALSE
  }
  if(orthogonalGP) {
    warning("X is rescaled to [-1,1]^p when orthogonalGP = TRUE\n")
    min.x <- apply(X, 2, min)
    range.x <- apply(X, 2, function(x) diff(range(x)))
    X <- t((t(X) - min.x)/range.x)
    X <- (X - 0.5) * 2
  }

  if(orthogonalGP & constantMean) {
    warning("orthogonalGP = FALSE since constantMean = TRUE\n")
    orthogonalGP <- FALSE
  }

  d <- ncol(X)
  if(is.null(d)){
    X <- matrix(X, ncol = 1)
    d <- 1
  }
  n <- nrow(X)
  N <- n * T_val
  #scale.x <- apply(X, 2, function(x) sqrt(sum(x^2)))
  scale.x <- apply(X, 2, function(x) max(abs(x)))
  X.scale <- t(t(X)/scale.x)    # scale X for theta_hat

  ###################       Model matrix       ###################
  y <- as.vector(Y)
  x.mx <- matrix(rep(t(X), T_val), ncol = d, byrow = TRUE)
  if(standardize) {
    mean.x <- apply(X, 2, mean)
    range.x <- apply(X, 2, function(x) diff(range(x)))
    x.mx <- t((t(x.mx) - mean.x)/range.x)
  }else{
    mean.x <- NULL
    range.x <- NULL
  }

  if(constantMean) M.mx <- matrix(1, ncol = 1, nrow = n) else M.mx <- cbind(1, x.mx)

  if(R > 0 | L > 0){
    y_lag.mx <- matrix(0, nrow = N, ncol = max(R,L))
    Ytemp <- Y
    for(i in 1:max(R,L)){
      Ytemp <- cbind(0, Ytemp)[,-(1+ncol(Ytemp))]
      y_lag.mx[,i] <- as.vector(Ytemp)
    }
  }
  if(R > 0){                          ### model matrix for y_lag
    M.mx <- cbind(M.mx, y_lag.mx[,1:R])
  }
  if(L > 0){                          ### model matrix for x, y interaction
    xy.mx <- matrix(0, nrow = N, ncol = d * L)
    for(i in 1:L) xy.mx[, (d *(i-1) + 1):(d * i)] <- x.mx * y_lag.mx[,i]
    #xy.mx <- xy.mx[,c(1,ncol(xy.mx)-1)]# temp
    M.mx <- cbind(M.mx, xy.mx)
  }

  ###################       Initialization       ###################

  sigma_hat <- 3                                              # initialize sigma_hat
  rho_hat <- rep(5, d)                                        # initialize theta_hat
  Z_hat <- rep(0, N)                                          # initialize Z_hat
  if(constantMean){                                           # initialize beta_hat
    fit0 <- glm(factor(y) ~ 1, family = binomial())
  }else{
    fit0 <- glm(factor(y) ~ M.mx - 1, family = binomial())
  }
  eta <- fit0$linear.predictors                               # initialize eta_hat
  eta_tilde <- eta + fit0$residuals                           # initialize eta_tilde_hat
  w <- fit0$weights                                           # initialize W

  ###################       Estimation (PQL/PQPL and REML)     ###################
  for (j in seq_len(maxit)) {

    eta.old <- eta
    ###################      Solve PQL/PQPL estimator     ###################

    if (verbose)
      message(gettextf("iteration %d", j), domain = NA)

    for (i in seq_len(maxit.PQPL)) {
      etaold <- eta
      R.mx <- corr_matrix(X, rho_hat, corr) + diag(nugget/(sigma_hat)^2, n)
      if(orthogonalGP) R.mx <- R.mx - orthogonalize(X, rho_hat, corr)
      out <- PQPL_estimate(M.mx, X, w, eta_tilde, R.mx, sigma_hat^2, T_val)
      eta <- out$eta
      beta_hat <- out$beta
      if (sum((eta - etaold)^2) < 1e-06 * sum(eta^2))
        break

      ### update ###
      p <- binomial()$linkinv(eta)
      p.var <- binomial()$variance(p)
      eta_tilde <- eta + (y - p)/p.var
      w <- p.var
    }

    if (sum((eta - eta.old)^2) < converge.tol * sum(eta^2))
      break

    ###################      Solve REML estimator     ###################
    if(is.null(nugget)){ ## with nugget term
      upper <- c(3, rep(3, d), 3)
      lower <- c(0.01, rep(-1.5, d), 0.01)

      if (j == 1){
        numIni <- 500
        para.int.mx <- t(t(maximinLHS(n = numIni, k = d + 2)) * (upper - lower) + lower)
        int.val <- rep(0, numIni)
        for(i in seq(numIni)) int.val[i] <- likelihood_w_nugget(para.int.mx[i,], M.mx = M.mx, X = X.scale, w = w, eta_tilde = eta_tilde, T_val = T_val, corr = corr, orthogonalGP = orthogonalGP)
        para.int.vt <- para.int.mx[which.min(int.val),]
      }

      para.vt <- lbfgs(para.int.vt, fn = likelihood_w_nugget, M.mx = M.mx, X = X.scale, w = w, eta_tilde = eta_tilde, T_val = T_val, corr = corr, orthogonalGP = orthogonalGP,
                       lower = lower, upper = upper, control = list(maxeval = maxit.REML, xtol_rel = xtol_rel))$par

      sigma_hat <- para.vt[1]
      if(corr$type == "exponential"){
        rho_hat <- para.vt[2:(d+1)] * log10(scale.x^corr$power)
      }else{
        rho_hat <- para.vt[2:(d+1)] * log10(scale.x)
      }
      nugget_hat <- para.vt[d+2]
      para.int.vt <- para.vt
    }else{
      upper <- c(3, rep(3, d))
      lower <- c(0.01, rep(-1.5, d))

      if (j == 1){
        numIni <- 500
        para.int.mx <- t(t(maximinLHS(n = numIni, k = d + 1)) * (upper - lower) + lower)
        int.val <- rep(0, numIni)
        for(i in seq(numIni)) int.val[i] <- likelihood(para.int.mx[i,], nugget = nugget, M.mx = M.mx, X = X.scale, w = w, eta_tilde = eta_tilde, T_val = T_val, corr = corr, orthogonalGP = orthogonalGP)
        para.int.vt <- para.int.mx[which.min(int.val),]
      }

      para.vt <- lbfgs(para.int.vt, fn = likelihood, nugget = nugget, M.mx = M.mx, X = X.scale, w = w, eta_tilde = eta_tilde, T_val = T_val, corr = corr, orthogonalGP = orthogonalGP,
                       lower = lower, upper = upper, control = list(maxeval = maxit.REML, xtol_rel = xtol_rel))$par

      sigma_hat <- para.vt[1]
      if(corr$type == "exponential"){
        rho_hat <- para.vt[2:(d+1)] + log10(scale.x^corr$power)
      }else{
        rho_hat <- para.vt[2:(d+1)] + log10(scale.x)
      }
      para.int.vt <- para.vt
    }
  }

  if(constantMean){
    alpha_hat <- beta_hat[1]
    if(R == 0) varphi_hat <- NULL else varphi_hat <- beta_hat[2:(1+R)]
    if(L == 0)  gamma_hat <- NULL else  gamma_hat <- beta_hat[(1+R):length(beta_hat)]
  }else{
    alpha_hat <- beta_hat[1:(d+1)]
    if(R == 0) varphi_hat <- NULL else varphi_hat <- beta_hat[(d+2):(d+1+R)]
    if(L == 0)  gamma_hat <- NULL else  gamma_hat <- beta_hat[(d+2+R):length(beta_hat)]
  }

  theta_hat <- 10^(-rho_hat)

  binaryGP <- list(alpha_hat = alpha_hat,
                   varphi_hat = varphi_hat,
                   gamma_hat = gamma_hat,
                   theta_hat = theta_hat,
                   sigma_hat = sigma_hat,
                   nugget_hat = ifelse(is.null(nugget), nugget_hat, nugget),
                   orthogonalGP = orthogonalGP,
                   X = X, Y = Y, R = R, L = L, corr = corr, constantMean = constantMean,
                   Model.mat = M.mx, eta_hat = eta,
                   standardize = standardize,
                   mean.x = mean.x,
                   range.x = range.x)

  class(binaryGP) <- "binaryGP"
  return(binaryGP)
}


