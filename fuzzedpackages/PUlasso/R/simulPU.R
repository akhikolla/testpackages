#' simulated PU data
#'
#' A simulated data for the illustration. 
#' Covariates \eqn{x_i} are drawn from \eqn{N(\mu,I_{5\times 5})}{N(\mu,I_{5x5})} or \eqn{N(-\mu,I_{5\times5})}{N(-\mu,I_{5x5})} with probability 0.5.
#' To make the first two variables active,\eqn{\mu = [\mu_1,\dots,\mu_2,0,0,0]^T, \theta = [\theta_0,\dots,\theta_2,0,0,0]^T}
#'  and we set \eqn{\mu_i=1.5, \theta_i \sim Unif[0.5,1]}{\mu_i=1.5, \theta_i ~ Unif[0.5,1]}
#' Responses \eqn{y_i} is simulated via \eqn{P_\theta(y=1|x) = 1/exp(-\theta^Tx)}.
#' 1000 observations are sampled from the sub-population of positives(y=1) and labeled, and another 1000 observations are sampled from the original population and unlabeled.
#'
#' @format A list containing model matrix X, true response y, labeled/unlabeled response vector z, and a true positive probability truePY1. 
"simulPU"
