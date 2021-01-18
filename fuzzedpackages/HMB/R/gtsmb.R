#' Generalized Two-Staged Model-Based estmation
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
#' @param Omega_S The covariance structure of \eqn{\boldsymbol{\epsilon}_{S}}{
#' \epsilon_S}, up to a constant.
#' @param Phis_Sa A 3D array, where the third dimension corresponds to the
#' covariance structure of
#' \eqn{E(\boldsymbol{\xi}_{k,Sa} \boldsymbol{\xi}_{j,Sa}^T)}{
#' E(\xi_k,Sa \xi_j,Sa')},
#' in the order \eqn{k=1, \ldots, p, j=1, \ldots k}{k=1,...,p, j=1,...,k}.
#' For p = 3, the order (k,j) will thus be (1,1), (2,1), (2,2), (3,1), (3,2),
#' (3,3).
#' @details
#' The GTSMB assumes the superpopulations
#' \deqn{y = \boldsymbol{x} \boldsymbol{\beta} + \epsilon}{
#'       y = x \beta + \epsilon}
#' \deqn{x_k = \boldsymbol{z} \boldsymbol{\gamma}_k + \xi_k}{
#'       x_k = z \gamma_k + \xi_k}
#' \deqn{\epsilon \perp \xi_k}{\epsilon indep. \xi_k}
#' For a sample from the superpopulation, the GTSMB assumes
#' \deqn{E(\boldsymbol{\epsilon}) = \mathbf{0},
#'       E(\boldsymbol{\epsilon} \boldsymbol{\epsilon}^T) = \omega^2 \boldsymbol{\Omega}}{
#'       E(\epsilon) = 0, E(\epsilon \epsilon') = \omega^2 \Omega}
#' \deqn{E(\boldsymbol{\xi}_k) = \mathbf{0},
#'       E(\boldsymbol{\xi}_k \boldsymbol{\xi}_j^T) = \theta_{\Phi,k,j}^2 \boldsymbol{\Phi}_{k,j},
#'       \theta_{\Phi,k,j}^2 \boldsymbol{\Phi}_{k,j} = \theta_{\Phi,j,k}^2 \boldsymbol{\Phi}_{j,k}}{
#'       E(\xi_k) = 0, E(\xi_k \xi_j') = \theta_Phi,k,j \Phi_k,j, \theta_Phi,k,j \Phi_k,j = \theta_Phi,j,k \Phi_j,k}
#' @return A fitted object of class HMB.
#' @seealso \code{\link{summary}},
#' \code{\link{getSpec}}.
#' @export
#' @examples
#' pop_U   = sample(nrow(HMB_data), 20000)
#' pop_Sa  = sample(pop_U, 500)
#' pop_S   = sample(pop_U, 100)
#'
#' y_S     = HMB_data[pop_S, "GSV"]
#' X_S     = HMB_data[pop_S, c("hMAX", "h80", "CRR")]
#' X_Sa    = HMB_data[pop_Sa, c("hMAX", "h80", "CRR")]
#' Z_Sa    = HMB_data[pop_Sa, c("B20", "B30", "B50")]
#' Z_U     = HMB_data[pop_U, c("B20", "B30", "B50")]
#'
#' Omega_S = diag(1, nrow(X_S))
#' Phis_Sa = array(0, c(nrow(X_Sa), nrow(X_Sa), ncol(X_Sa) * (ncol(X_Sa) + 1) / 2))
#' Phis_Sa[, , 1] = diag(1, nrow(X_Sa)) # Phi(1,1)
#' Phis_Sa[, , 2] = diag(1, nrow(X_Sa)) # Phi(2,1)
#' Phis_Sa[, , 3] = diag(1, nrow(X_Sa)) # Phi(2,2)
#' Phis_Sa[, , 4] = diag(1, nrow(X_Sa)) # Phi(3,1)
#' Phis_Sa[, , 5] = diag(1, nrow(X_Sa)) # Phi(3,2)
#' Phis_Sa[, , 6] = diag(1, nrow(X_Sa)) # Phi(3,3)
#'
#' gtsmb_model = gtsmb(y_S, X_S, X_Sa, Z_Sa, Z_U, Omega_S, Phis_Sa)
#' gtsmb_model
#' @references Holm, S., Nelson, R. & Ståhl, G. (2017) Hybrid three-phase estimators for large-area forest inventory using ground plots, 
#' airborne lidar, and space lidar. \emph{Remote Sensing of Environment, 197,} 85–97.
#' 
#' Saarela, S., Holm, S., Healey, S.P., Andersen, H.-E., Petersson, H., Prentius, W., Patterson, P.L., Næsset, E., Gregoire, T.G. & Ståhl, G. (2018). 
#' Generalized Hierarchical Model-Based Estimation for Aboveground Biomass Assessment Using GEDI and Landsat Data, \emph{Remote Sensing, 10(11),} 1832.

#' @export
gtsmb = function(
  y_S,
  X_S,
  X_Sa,
  Z_Sa,
  Z_U,
  Omega_S,
  Phis_Sa
) {
  popCheck(y_S, X_S, X_Sa, Z_Sa, Z_U)

  if (
    missingArg(Omega_S) ||
    missingArg(Phis_Sa)
  ) {
    stop('Missing necessary arguments.')
  }

  if (
    nrow(X_S) != nrow(Omega_S) ||
    nrow(Omega_S) != ncol(Omega_S)
  ) {
    stop('Omega_S has incorrect dimensions.')
  }

  if (
    !is.array(Phis_Sa) ||
    dim(Phis_Sa)[3] != ncol(X_Sa) * (ncol(X_Sa) + 1) / 2
  ) {
    stop('Phis_Sa is of incorrect type or size.')
  }

  ## Initialize HMB
  h = new("HMB")
  h@method = 'GTSMB'

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
    Phis_Sa = as.array(Phis_Sa)
  )

  model = cpp_gtsmb(
    h@data$y_S,
    h@data$X_S,
    h@data$X_Sa,
    h@data$Z_Sa,
    h@data$Z_U,
    h@data$Omega_S,
    h@data$Phis_Sa
  )

  h@Alpha = model$Gamma %*% model$Beta
  h@AlphaCov = model$GammaCov_ish
  h@Gamma = model$Gamma
  h@Beta = model$Beta
  h@BetaCov = model$BetaCov
  h@mu = model$mu
  h@muVar = model$muVar
  
  h@resids = list()
    h@resids$omega2 = model$omega
    h@resids$phi2s = model$phi2s
  

  return(h)
}
