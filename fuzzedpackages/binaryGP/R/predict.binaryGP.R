#' Predictions of Binary Gaussian Process
#'
#' @description The function computes the predicted response and its variance as well as its confidence interval.
#'
#' @param object a class binaryGP object estimated by \code{binaryGP_fit}.
#' @param xnew a testing matrix with dimension \code{n_new} by \code{d} in which each row corresponds to a predictive location.
#' @param conf.level a value from 0 to 1 specifying the level of confidence interval. The default is 0.95.
#' @param sim.number a positive integer specifying the simulation number for Monte-Carlo method. The default is 101.
#' @param ... for compatibility with generic method \code{predict}.
#'
#' @return
#' \item{mean}{a matrix with dimension \code{n_new} by \code{T} displaying predicted responses at locations \code{xnew}.}
#' \item{var}{a matrix with dimension \code{n_new} by \code{T} displaying predictive variances at locations \code{xnew}.}
#' \item{upper.bound}{a matrix with dimension \code{n_new} by \code{T} displaying upper bounds with \code{conf.level} confidence level.}
#' \item{lower.bound}{a matrix with dimension \code{n_new} by \code{T} displaying lower bounds with \code{conf.level} confidence level.}
#' \item{y_pred}{a matrix with dimension \code{n_new} by \code{T} displaying predicted binary responses at locations \code{xnew}.}
#'
#' @seealso \code{\link{binaryGP_fit}} for estimation of the binary Gaussian process.
#' @author Chih-Li Sung <iamdfchile@gmail.com>
#'
#' @import logitnorm
#' @examples
#'
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
#' X.test <- matrix(runif(d * n.test), ncol = d)
#'
#' ##### without time-series #####
#' Y <- test_function(X, 1)$Y  ## Y is a vector
#' test.out <- test_function(X.test, 1)
#' Y.test <- test.out$Y
#' P.true <- test.out$P
#' \donttest{
#' # fitting
#' binaryGP.model <- binaryGP_fit(X = X, Y = Y)
#'
#' # prediction
#' binaryGP.prediction <- predict(binaryGP.model, xnew = X.test)
#' print(binaryGP.prediction$mean)
#' print(binaryGP.prediction$var)
#' print(binaryGP.prediction$upper.bound)
#' print(binaryGP.prediction$lower.bound)
#'
#' ##### with time-series #####
#' Y <- test_function(X, 10)$Y  ## Y is a matrix with 10 columns
#' test.out <- test_function(X.test, 10)
#' Y.test <- test.out$Y
#' P.true <- test.out$P
#'
#' # fitting
#' binaryGP.model <- binaryGP_fit(X = X, Y = Y, R = 1)
#'
#' # prediction
#' binaryGP.prediction <- predict(binaryGP.model, xnew = X.test)
#' print(binaryGP.prediction$mean)
#' print(binaryGP.prediction$var)
#' print(binaryGP.prediction$upper.bound)
#' print(binaryGP.prediction$lower.bound)
#' }
#' @export predict.binaryGP
#' @export
#'

predict.binaryGP <- function(object, xnew, conf.level = 0.95, sim.number = 101, ...){

  if (is.binaryGP(object) == FALSE) {
    stop("The object in question is not of class binaryGP \n")
  }
  if (is.matrix(xnew) == FALSE) {
    xnew <- as.matrix(xnew)
  }

  X            <- object$X
  Y            <- object$Y
  R            <- object$R
  L            <- object$L
  corr         <- object$corr
  constantMean <- object$constantMean
  alpha_hat    <- object$alpha_hat
  varphi_hat   <- object$varphi_hat
  gamma_hat    <- object$gamma_hat
  theta_hat    <- object$theta_hat
  rho_hat      <- -log10(theta_hat)
  sigma_hat    <- object$sigma_hat
  nugget       <- object$nugget_hat
  eta_hat      <- object$eta_hat
  M.mx         <- object$Model.mat
  standardize  <- object$standardize
  orthogonalGP <- object$orthogonalGP
  mean.x       <- object$mean.x
  range.x      <- object$range.x

  ###################       Setting       ###################
  n <- nrow(X)
  n.test <- nrow(xnew)
  d <- ncol(X)
  T_val <- ncol(Y)
  if(is.null(T_val)) {
    T_val <- 1
    Y <- matrix(Y, ncol = 1)
  }
  N <- n * T_val

  if(standardize){
    xnew <- t((t(xnew) - mean.x)/range.x)
  }

  R.mx <- corr_matrix(X, rho_hat, corr)
  R.mx <- R.mx + diag(nugget/(sigma_hat)^2, n)  # add a nugget term to make the variance-covariance matrix more stable
  if(orthogonalGP) R.mx <- R.mx - orthogonalize(X, rho_hat, corr)
  R_inv <- solve(R.mx)

  p.vt <- binomial()$linkinv(eta_hat)
  p.vt[p.vt > 1-1e-8] <- 1-1e-8
  p.vt[p.vt < 1e-8] <- 1e-8
  p.mx <- matrix(p.vt, ncol = T_val)  ## initial p

  beta_hat <- c(alpha_hat, varphi_hat, gamma_hat)
  mu <- c(M.mx %*% beta_hat)
  mu.mx <- matrix(mu, ncol = T_val)

  P.new <- P_mean.new <- P_var.new <-
    P_quantile_u.new <- P_quantile_l.new <- Y.new <- array(NA, c(n.test, sim.number, T_val))

  for(j in 1:sim.number){
    ###################      MH algorithm     ###################
    for(tt in 1:T_val){
      for(k in 1:n){
        Q_kk <- R_inv[k,k]
        Q_kj <- R_inv[k,-k]
        m_D <- mu.mx[k,tt] - (sum(Q_kj * (log(p.mx[-k,tt]) - log(1-p.mx[-k,tt]) - mu.mx[-k,tt])))/Q_kk
        v_D <- sigma_hat^2/Q_kk

        p_star <- rlogitnorm(n = 1, mu = m_D, sigma = sqrt(v_D))
        p_star <- min(c(p_star, 1-1e-8))
        p_star <- max(c(p_star, 1e-8))
        U <- runif(1, 0, 1)
        f_1 <- dbinom(Y[k,tt], size = 1, prob = p_star)
        f_2 <- dbinom(Y[k,tt], size = 1, prob = p.mx[k,tt])
        if(U < min(c(1, f_1/f_2))){
          p.mx[k,tt] <- p_star
        }
      }
    }
    ###################      w/o time-series: monte carlo method                                           ###################
    ###################      w/ time-series : draw conditional distribution and then monte carlo method    ###################
    for(i in 1:n.test){
      x.new <- xnew[i,,drop=FALSE]
      p.new <- p_mean.new <- p_var.new <-
        p_quantile_u.new  <- p_quantile_l.new <- y.new <- rep(NA, T_val)
      r.vt <- corr_vec(x.new, X, rho_hat, corr)
      if(orthogonalGP) r.vt <- r.vt - orthogonalize_vec(x.new, X, rho_hat, corr)
      k.vt <- r.vt %*% R_inv
      rRinvr <- c(k.vt %*% t(r.vt))
      for(tt in 1:T_val){
        if(constantMean)  M.mx.new <- matrix(1, ncol = 1, nrow = 1) else M.mx.new <- cbind(1, x.new)
        if(R > 0){
          if(tt <= R){
            ytmp <- rep(0,R)
            ytmp[0:length(y.new[(tt-1):max(tt-R, 0)])] <- y.new[(tt-1):max(tt-R, 0)]
          }else{
            ytmp <- y.new[(tt-1):(tt-R)]
          }
          M.mx.new <- cbind(M.mx.new, matrix(ytmp, nrow=1))
        }
        if(L > 0){
          if(tt <= L){
            ytmp <- rep(0,L)
            ytmp[0:length(y.new[(tt-1):max(tt-L, 0)])] <- y.new[(tt-1):max(tt-L, 0)]
          }else{
            ytmp <- y.new[(tt-1):(tt-L)]
          }
          xy.tmp <- matrix(0, nrow = 1, ncol = d * L)
          for(l in 1:L) xy.tmp[1,(d *(l-1) + 1):(d * l)] <- x.new * ytmp[l]
          #xy.tmp <- xy.tmp[,c(ncol(xy.tmp)-1), drop=FALSE] #temp
          M.mx.new <- cbind(M.mx.new, xy.tmp)
        }

        m_D <- c(M.mx.new %*% beta_hat) + k.vt %*% (log(p.mx[,tt])-log(1-p.mx[,tt])-mu.mx[,tt])
        v_D <- sigma_hat^2 * max((1 - rRinvr), 0)

        p_mean.new[tt] <- momentsLogitnorm(mu = m_D, sigma = sqrt(v_D), stop.on.error = FALSE)[1]
        p_var.new[tt] <- momentsLogitnorm(mu = m_D, sigma = sqrt(v_D), stop.on.error = FALSE)[2]

        if(sqrt(v_D) < 1e-3 | p_mean.new[tt] < 1e-4){ ### bug from momentsLogitnorm package due to sigma is too small
          p_mean.new[tt] <- rlogitnorm(n = 1, mu = m_D, sigma = sqrt(v_D))
        }

        p_quantile_u.new[tt] <- qlogitnorm(p = 1 - (1 - conf.level)/2, mu = m_D, sigma = sqrt(v_D))
        p_quantile_l.new[tt] <- qlogitnorm(p = (1 - conf.level)/2, mu = m_D, sigma = sqrt(v_D))

        if(T_val > 1){
          ### update if time-series is included ###
          p.new[tt] <- rlogitnorm(n = 1, mu = m_D, sigma = sqrt(v_D))
          y.new[tt] <- rbinom(n = 1, size = 1, prob = p.new[tt])
        }else{
          p.new[tt] <- rlogitnorm(n = 1, mu = m_D, sigma = sqrt(v_D))
          y.new[tt] <- rbinom(n = 1, size = 1, prob = p_mean.new[tt])
        }
      }
      P_mean.new[i, j, ] <- p_mean.new
      P_var.new[i, j, ] <- p_var.new
      P_quantile_u.new[i, j, ] <- p_quantile_u.new
      P_quantile_l.new[i, j, ] <- p_quantile_l.new
      P.new[i, j, ] <- p.new
      Y.new[i, j, ] <- y.new
    }
  }


  if(T_val > 1){
    p_new <- apply(P.new, c(1,3), median)
    p_var_new <- apply(P.new, c(1,3), var)
    p_quantile_u_new <- apply(P.new, c(1,3), quantile, probs = 1 - (1 - conf.level)/2)
    p_quantile_l_new <- apply(P.new, c(1,3), quantile, probs = (1 - conf.level)/2)
  }else{
    p_new <- apply(P_mean.new, c(1,3), median)
    p_var_new <- apply(P_var.new, c(1,3), median) + apply(P_mean.new, c(1,3), var)
    p_quantile_u_new <- apply(P.new, c(1,3), quantile, probs = 1 - (1 - conf.level)/2)
    p_quantile_l_new <- apply(P.new, c(1,3), quantile, probs = (1 - conf.level)/2)
  }

  y_new <- apply(Y.new, c(1,3), median)

  return(list(mean = p_new, var = p_var_new, upper.bound = p_quantile_u_new, lower.bound = p_quantile_l_new, y_pred = y_new))
}
