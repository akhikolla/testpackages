#' @title Calculate Theoretical Allan Variance for Stationary First-Order Autoregressive
#' (AR1) Process
#' @description
#' This function allows us to calculate the theoretical allan variance for stationary
#' first-order autoregressive (AR1) process.
#' @export
#' @param n An \code{integer} value for the size of the cluster.
#' @param phi A \code{double} value for the autocorrection parameter \eqn{\phi}{phi}.
#' @param sigma2 A \code{double} value for the variance parameter \eqn{\sigma ^2}{sigma^2}.
#' @return A \code{double} indicating the theoretical allan variance for AR1 process.
#' @note This function is based on the calculation of the theoretical allan variance
#' for stationary AR1 process raised in "Allan Variance of Time Series Models for
#' Measurement Data" by Nien Fan Zhang, 2008, Metrologia, 45(5): 549.
#' This calculation is fundamental and necessary for the study in "A Study of the Allan Variance
#' for Constant-Mean Non-Stationary Processes" by Xu et al., 2017, IEEE Signal Processing Letters, 24(8): 1257–1260.
#' @author Yuming Zhang
#' @examples
#' av1 = av_ar1(n = 5, phi = 0.9, sigma2 = 1)
#' av2 = av_ar1(n = 8, phi = 0.5, sigma2 = 2)
av_ar1 = function(n, phi, sigma2){
  numerator = n-3*phi-n*phi^2+4*phi^(n+1)-phi^(2*n+1)
  denominator = n^2 * (1-phi)^2 * (1-phi^2)
  result = numerator / denominator * sigma2
  return(result)
}




#' @title Calculate Theoretical Allan Variance for Stationary White Noise Process
#' @description
#' This function allows us to calculate the theoretical allan variance for stationary
#' white noise process.
#' @export
#' @param sigma2 A \code{double} value for the variance parameter \eqn{\sigma ^2}{sigma^2}.
#' @param n An \code{integer} value for the size of the cluster.
#' @return A \code{double} indicating the theoretical allan variance for the white noise
#' process.
#' @note This function is based on the calculation of the theoretical allan variance
#' for stationary white noise process raised in "Allan Variance of Time Series Models for
#' Measurement Data" by Nien Fan Zhang, 2008, Metrologia, 45(5): 549.
#' This calculation is fundamental and necessary for the study in "A Study of the Allan Variance
#' for Constant-Mean Non-Stationary Processes" by Xu et al., 2017, IEEE Signal Processing Letters, 24(8): 1257–1260.
#' @examples
#' av1 = av_wn(sigma2 = 1, n = 5)
#' av2 = av_wn(sigma2 = 2, n = 8)
av_wn = function(sigma2, n){
  result = sigma2/n
  return(result)
}


#' @title Calculate Theoretical Allan Variance for Stationary Quantization Noise Process
#' @description
#' This function allows us to calculate the theoretical allan variance for stationary
#' quantization noise process.
#' @export
#' @param Q2 A \code{double} value for the noise parameter \eqn{Q^2}{Q^2}.
#' @param n An \code{integer} value for the size of the cluster.
#' @return A \code{double} indicating the theoretical allan variance for the quantization noise
#' process.
#' @note This function is based on the calculation of the theoretical allan variance
#' for stationary quantization noise process raised in "Allan Variance of Time Series Models for
#' Measurement Data" by Nien Fan Zhang, 2008, Metrologia, 45(5): 549.
#' This calculation is fundamental and necessary for the study in "A Study of the Allan Variance
#' for Constant-Mean Non-Stationary Processes" by Xu et al., 2017, IEEE Signal Processing Letters, 24(8): 1257–1260.
#' @examples
#' av1 = av_qn(Q2 = 1, n = 5)
#' av2 = av_qn(Q2 = 2, n = 8)
av_qn = function(Q2, n){
  result = 3*Q2 / n^2
  return(result)
}


#' @title Calculate Theoretical Allan Variance for Random Walk Process
#' @description
#' This function allows us to calculate the theoretical allan variance for
#' random walk process.
#' @export
#' @param omega2 A \code{double} value for the noise parameter \eqn{\omega ^2}{omega^2}.
#' @param n An \code{integer} value for the size of the cluster.
#' @return A \code{double} indicating the theoretical allan variance for the random walk
#' process.
#' @note This function is based on the calculation of the theoretical allan variance
#' for random walk process raised in "Allan Variance of Time Series Models for
#' Measurement Data" by Nien Fan Zhang, 2008, Metrologia, 45(5): 549.
#' This calculation is fundamental and necessary for the study in "A Study of the Allan Variance
#' for Constant-Mean Non-Stationary Processes" by Xu et al., 2017, IEEE Signal Processing Letters, 24(8): 1257–1260.
#' @examples
#' av1 = av_rw(omega2 = 1, n = 5)
#' av2 = av_rw(omega2 = 2, n = 8)
av_rw = function(omega2, n){
  result = (2 * n^2 + 1)*omega2 / (6*n)
  return(result)
}


#' @title Calculate Theoretical Allan Variance for Drift Process
#' @description
#' This function allows us to calculate the theoretical allan variance for
#' drift process.
#' @export
#' @param delta A \code{double} value for the noise parameter \eqn{\delta}{delta}.
#' @param n An \code{integer} value for the size of the cluster.
#' @return A \code{double} indicating the theoretical allan variance for the drift
#' process.
#' @note This function is based on the calculation of the theoretical allan variance
#' for drift process raised in "Allan Variance of Time Series Models for
#' Measurement Data" by Nien Fan Zhang, 2008, Metrologia, 45(5): 549.
#' This calculation is fundamental and necessary for the study in "A Study of the Allan Variance
#' for Constant-Mean Non-Stationary Processes" by Xu et al., 2017, IEEE Signal Processing Letters, 24(8): 1257–1260.
#' @examples
#' av1 = av_dr(delta = 1, n = 5)
#' av2 = av_dr(delta = 2, n = 8)
av_dr = function(delta, n){
  result = (n^2 * delta^2) / 2
  return(result)
}
