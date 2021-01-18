#' Generalized Hierarchical Model-Based estimation method
#' @encoding UTF-8
#' @param y_S Response object that can be coerced into a column vector. The
#' \code{_S} denotes that \code{y} is part of the sample \emph{S}, with
#' \eqn{N_S \le N_{Sa} \le N_U}{N_S \le N_Sa \le N_U}.
#' @param X_S Object of predictors variables that can be coerced into a matrix.
#' The rows of \code{X_S} correspond to the rows of \code{y_S}.
#' @param X_Sa Object of predictor variables that can be coerced into a matrix.
#' The set \emph{Sa} is the intermediate sample.
#' @param Z_Sa Object of predictor variables that can be coerced into a matrix.
#' The set \emph{Sa} is the intermediate sample, and the Z-variables often some
#' sort of auxilairy, inexpensive data. The rows of \code{Z_Sa} correspond to
#' the rows of \code{X_Sa}
#' @param Z_U Object of predictor variables that can be coerced into a matrix.
#' The set \emph{U} is the universal population sample.
#' @param Omega_S The covariance structure of \eqn{\boldsymbol{\epsilon}_{S}}{
#' \epsilon_S}, up to a constant.
#' @param Sigma_Sa The covariance structure of \eqn{\boldsymbol{u}_{Sa}}{
#' u_Sa}, up to a constant.
#' @details
#' The GHMB assumes two models
#' \deqn{y = \boldsymbol{x} \boldsymbol{\beta} + \epsilon}{
#'       y = x \beta + \epsilon}
#' \deqn{\boldsymbol{x} \boldsymbol{\beta} = \boldsymbol{z} \boldsymbol{\alpha} + \boldsymbol{u}}{
#'       x \beta = z \alpha + u}
#' \deqn{\epsilon \perp u}{\epsilon indep. u}
#' For a sample from the superpopulation, the GHMB assumes
#' \deqn{E(\boldsymbol{\epsilon}) = \mathbf{0},
#'       E(\boldsymbol{\epsilon} \boldsymbol{\epsilon}^T) = \omega^2 \boldsymbol{\Omega}}{
#'       E(\epsilon) = 0, E(\epsilon \epsilon') = \omega^2 \Omega}
#' \deqn{E(\boldsymbol{u}) = \mathbf{0},
#'       E(\boldsymbol{u} \boldsymbol{u}^T) = \sigma^2 \boldsymbol{\Sigma}}{
#'       E(u) = 0, E(u u') = \sigma^2 \Sigma}
#' @return A fitted object of class HMB. 
#' @seealso \code{\link{summary}},
#' \code{\link{getSpec}}.
#' @export
#' @examples
#' pop_U    = sample(nrow(HMB_data), 20000)
#' pop_Sa   = sample(pop_U, 2500)
#' pop_S    = sample(pop_U, 300)
#'
#' y_S      = HMB_data[pop_S, "GSV"]
#' X_S      = HMB_data[pop_S, c("hMAX", "h80", "CRR", "pVeg")]
#' X_Sa     = HMB_data[pop_Sa, c("hMAX", "h80", "CRR", "pVeg")]
#' Z_Sa     = HMB_data[pop_Sa, c("B20", "B30", "B50")]
#' Z_U      = HMB_data[pop_U, c("B20", "B30", "B50")]
#'
#' Omega_S  = diag(1, nrow(X_S))
#' Sigma_Sa = diag(1, nrow(Z_Sa))
#'
#' ghmb_model = ghmb(
#'   y_S, X_S, X_Sa, Z_Sa, Z_U, Omega_S, Sigma_Sa)
#' ghmb_model
#' @references Saarela, S., Holm, S., Healey, S.P., Andersen, H.-E., Petersson, H., Prentius, W., Patterson, P.L., Næsset, E., Gregoire, T.G. & Ståhl, G. (2018). 
#' Generalized Hierarchical Model-Based Estimation for Aboveground Biomass Assessment Using GEDI and Landsat Data, \emph{Remote Sensing, 10(11),} 1832.

#' @export
ghmb = function(
  y_S,
  X_S,
  X_Sa,
  Z_Sa,
  Z_U,
  Omega_S,
  Sigma_Sa) {
  popCheck(y_S, X_S, X_Sa, Z_Sa, Z_U)

  if (
    missingArg(Omega_S) ||
    missingArg(Sigma_Sa)
  ) {
    stop('Missing necessary arguments.')
  }

  if (
    nrow(Z_Sa) != nrow(Sigma_Sa) ||
    nrow(Sigma_Sa) != ncol(Sigma_Sa)
  ) {
    stop('Sigma_Sa has incorrect dimensions.')
  }

  if (
    nrow(X_S) != nrow(Omega_S) ||
    nrow(Omega_S) != ncol(Omega_S)
  ) {
    stop('Omega_S has incorrect dimensions.')
  }



  ## Initialize HMB
  h = new("HMB")
  h@method = 'GHMB'

  h@n = list(
    U = nrow(Z_U),
    Sa = nrow(X_Sa),
    S = nrow(X_S)
  )

  h@data = list(
    y_S = as.matrix(y_S),
    X_S = as.matrix(cbind(rep(1, h@n$S), X_S)),
    X_Sa = as.matrix(cbind(rep(1, h@n$Sa), X_Sa)),
    Z_Sa = as.matrix(cbind(rep(1, h@n$Sa), Z_Sa)),
    Z_U = as.matrix(cbind(rep(1, h@n$U), Z_U)),
    Omega_S = as.matrix(Omega_S),
    Sigma_Sa = as.matrix(Sigma_Sa)
  )

  model = cpp_ghmb(
    h@data$y_S,
    h@data$X_S,
    h@data$X_Sa,
    h@data$Z_Sa,
    h@data$Z_U,
    h@data$Omega_S,
    h@data$Sigma_Sa)

  h@resids = list(
    SigmaConst = model$sigma,
    OmegaConst = model$omega
  )

  h@Alpha = model$Alpha
  h@AlphaCov = model$AlphaCov
  h@Beta = model$Beta
  h@BetaCov = model$BetaCov
  h@mu = model$mu
  h@muVar = model$muVar

  return(h)
}
