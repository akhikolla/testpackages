#' Summary of Fitting a Binary Gaussian Process
#'
#' @description The function summarizes estimation and significance results by \code{binaryGP_fit}.
#'
#' @param object a class binaryGP object estimated by \code{binaryGP_fit}.
#' @param ... for compatibility with generic method \code{summary}.
#'
#' @return A table including the estimates by \code{binaryGP_fit}, and the correponding standard deviations, Z-values and p-values.
#' @seealso \code{\link{binaryGP_fit}} for estimation of the binary Gaussian process.
#' @author Chih-Li Sung <iamdfchile@gmail.com>
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
#' @export

summary.binaryGP <- function(object, ...){

  if (is.binaryGP(object) == FALSE) {
    stop("The object in question is not of class binaryGP \n")
  }
  if (object$constantMean){
    stop("Not available when object$constantMean = TRUE \n")
  }

  X           <- object$X
  Y           <- object$Y
  R           <- object$R
  L           <- object$L
  alpha_hat   <- object$alpha_hat
  varphi_hat  <- object$varphi_hat
  gamma_hat   <- object$gamma_hat
  theta_hat   <- object$theta_hat
  sigma_hat   <- object$sigma_hat
  eta_hat     <- object$eta_hat
  M.mx        <- object$Model.mat

  ###################       Setting       ###################
  n <- nrow(X)
  d <- ncol(X)
  T_val <- ncol(Y)
  if(is.null(T_val)) {
    T_val <- 1
    Y <- matrix(Y, ncol = 1)
  }
  N <- n * T_val
  beta_hat <- c(alpha_hat, varphi_hat, gamma_hat)

  ###################       Variance of estimators       ###################
  p.vt <- binomial()$linkinv(eta_hat)
  #Gamma_N <- t(M.mx) %*% diag(c(p.vt * (1 - p.vt))) %*% M.mx / N
  Gamma_N <- t(M.mx) %*% (M.mx*c(p.vt * (1 - p.vt))) / N
  Cov_beta <- solve(Gamma_N) / N
  var_beta <- diag(Cov_beta)

  ###################       Summary table      ###################
  col.names <- c("Intercept", paste0("alpha_", 1:(length(alpha_hat)-1)))
  if(R > 0) col.names <- c(col.names, paste0("varphi_", 1:length(varphi_hat)))
  if(L > 0) col.names <- c(col.names, paste0("gamma_", 1:length(gamma_hat)))
  row.names <- c("Coefficients", "Standard_deviations", "Z_values", "P_values")

  summary.tbl <- matrix(0, nrow = length(beta_hat), ncol = 4,
                        dimnames = list(col.names, row.names))
  summary.tbl[,"Coefficients"] <- beta_hat
  summary.tbl[,"Standard_deviations"] <- sqrt(var_beta)
  summary.tbl[,"Z_values"] <- beta_hat/sqrt(var_beta)
  summary.tbl[,"P_values"] <- (1 - pnorm(abs(beta_hat)/sqrt(var_beta))) * 2

  summary.df <- data.frame(summary.tbl)
  summary.df[,"P_values"] <- format(round(summary.df[,"P_values"],4), digits = 4)

  return(summary.df)
}
