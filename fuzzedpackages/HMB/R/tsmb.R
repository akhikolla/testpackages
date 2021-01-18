#' Two-staged Model-Based estmation
#' @encoding UTF-8
#' @param y_S Response object that can be coersed into a column vector. The
#' \code{_S} denotes that \code{y} is part of the sample \emph{S}, with
#' \eqn{N_S \le N_{Sa} \le N_U}{N_S \le N_Sa \le N_U}.
#' @param X_S Object of predictors variables that can be coersed into a matrix.
#' The rows of \code{X_S} correspond to the rows of \code{y_S}.
#' @param X_Sa Object of predictor variables that can be coresed into a matrix.
#' The set \emph{Sa} is the intermediate sample.
#' @param Z_Sa Object of predictor variables that can be coresed into a matrix.
#' The set \emph{Sa} is the intermediate sample, and the Z-variables often some
#' sort of auxilairy, inexpensive data. The rows of \code{Z_Sa} correspond to
#' the rows of \code{X_Sa}
#' @param Z_U Object of predictor variables that can be coresed into a matrix.
#' The set \emph{U} is the universal population sample.
#' @return A fitted object of class HMB.
#' @details
#' The TSMB assumes the superpopulations
#' \deqn{y = \boldsymbol{x}^T \boldsymbol{\beta} + \epsilon}{
#'       y = x' \beta + \epsilon}
#' \deqn{x_k = \boldsymbol{z}^T \boldsymbol{\gamma}_k + \xi_k}{
#'       x_k = z' \gamma_k + \xi_k}
#' \deqn{\epsilon \perp \xi_k}{\epsilon indep. \xi_k}
#' For a sample from the superpopulation, the TSMB assumes
#' \deqn{E(\boldsymbol{\epsilon}) = \mathbf{0},
#'       E(\boldsymbol{\epsilon} \boldsymbol{\epsilon}^T) = \omega^2 \mathbf{I}}{
#'       E(\epsilon) = 0, E(\epsilon \epsilon') = \omega^2 I}
#' \deqn{E(\boldsymbol{\xi}_k) = \mathbf{0},
#'       E(\boldsymbol{\xi}_k \boldsymbol{\xi}_j^T) = \phi_{k,j}^2 \mathbf{I}}{
#'       E(\xi_k) = 0, E(\xi_k \xi_j') = \phi_k,j^2 I}
#' @seealso \code{\link{summary}},
#' \code{\link{getSpec}}.
#' @export
#' @examples
#' pop_U  = sample(nrow(HMB_data), 20000)
#' pop_Sa = sample(pop_U, 5000)
#' pop_S  = sample(pop_U, 300)
#'
#' y_S    = HMB_data[pop_S, "GSV"]
#' X_S    = HMB_data[pop_S, c("hMAX", "h80", "CRR", "pVeg")]
#' X_Sa   = HMB_data[pop_Sa, c("hMAX", "h80", "CRR", "pVeg")]
#' Z_Sa   = HMB_data[pop_Sa, c("B20", "B30", "B50")]
#' Z_U    = HMB_data[pop_U, c("B20", "B30", "B50")]
#'
#' tsmb_model = tsmb(y_S, X_S, X_Sa, Z_Sa, Z_U)
#' tsmb_model
#' @references Saarela, S., Holm, S., Grafström, A., Schnell, S., Næsset, E., Gregoire, T.G., Nelson, R.F. & Ståhl, G. (2016). 
#' Hierarchical model-based inference for forest inventory utilizing three sources of information. \emph{Annals of Forest Science, 73(4),} 895-910.

#' @export
tsmb = function(
  y_S,
  X_S,
  X_Sa,
  Z_Sa,
  Z_U) {
  popCheck(y_S, X_S, X_Sa, Z_Sa, Z_U)

  h = new("HMB")
  h@method = 'TSMB'

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
    Z_U = as.matrix(cbind(rep(1, h@n$U), Z_U))
  )

  model = cpp_tsmb(
    h@data$y_S,
    h@data$X_S,
    h@data$X_Sa,
    h@data$Z_Sa,
    h@data$Z_U)

  h@resids = list()
    h@resids$omega2 = model$omega
    h@resids$phi2s = model$phis2
    h@resids$sigma2 = model$sigma2

  h@Alpha = model$Gamma %*% model$Beta
  h@AlphaCov = model$AlphaCov
  h@Beta = model$Beta
  h@BetaCov = model$BetaCov
  h@Gamma = model$Gamma

  h@mu = model$mu
  h@muVar = model$muVar

  return(h)
}
