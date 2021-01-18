#' Bayesian wavelet denoising
#'
#' This function carries out Bayesian wavelet denoising using the Normal 
#' Inverse Gamma Markov Tree method of Ma and Soriano (2016).
#'
#' @param W An object of class \code{DWT}.
#' @param alpha Hyperparameter controlling the global smoothness.
#' @param nu Hyperparameter controlling variance heterogeneity. If \code{Inf},
#' then the variance is identical for all nodes.
#' @param n.samples Number of posterior draws.
#' @param transition.mode Type of transition. 
#' The two options are \code{Markov} or \code{Independent}.
#' @param method Method used for find maxmimum of marginal likelihood. 
#' 
#' @return An object of class \code{grove}.
#' @references Ma L. and Soriano J. (2016) Efficient functional ANOVA 
#' through wavelet-domain Markov groves. arXiv:1602.03990v2 [stat.ME]
#' (\url{https://arxiv.org/abs/1602.03990v2}).
#' @export
#' @examples
#' data <- wavethresh::DJ.EX(n = 512, noisy = TRUE, rsnr = 5)$doppler
#' W <- DWT(data)
#' ans <- Denoise(W)
Denoise <- function(W, 
                    alpha = 0.5, 
                    nu = 5, 
                    n.samples = 500, 
                    transition.mode = "Markov",
                    method = "Nelder-Mead") {
  
  frml <- formula(~ 1)
  X <- data.frame(rep(1, nrow(W$D)))
  
  output <- .groveEB.all.random(W = W, 
                                formula = frml, 
                                X = X, 
                                alpha = alpha, 
                                nu = nu, 
                                eta.rho = .InitEtaRhoPar(),
                                eta.kappa = .InitEtaKappaPar(),
                                gamma.rho = .InitGammaRhoPar(),
                                gamma.kappa = .InitGammaKappaPar(),
                                n.samples = n.samples, 
                                verbose = FALSE,
                                transition.mode = transition.mode,
                                method = method)
  return(output)
}  
