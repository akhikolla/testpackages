#' @useDynLib gfiExtremes
#' @importFrom Rcpp evalCpp
#' @import coda
NULL

rGPDmodel <- function(nsims, mu = 10, gamma = 1, sigma = 1){
  U <- runif(nsims)
  P <- qgpareto(1-U, mu, gamma, sigma)
  p <- 1 / (1/mu + 1)   # p * 1/a = 1 - p => p = 1/(1/a+1)
  X <- numeric(nsims)
  for(i in 1L:nsims){
    X[i] <- sample(c(mu*U[i], P[i]), 1L, prob = c(p, 1-p))
  }
  X
}

thinChain <- function(chain, skip){
  niterations <- nrow(chain)
  every.ith <- c(TRUE, rep(FALSE, skip))
  keep <- rep(every.ith, ceiling(niterations / (skip+1L)))[1L:niterations]
  chain[keep, ]
}

thresholdIndex <- function(Xs, a){
  # Xs is X sorted
  match(TRUE, Xs - a >= 0)
}

#' @importFrom stats optim var
#' @noRd
gpdFit <- function(X, mu){
  fn <- function(gamma_sigma){
    gamma <- gamma_sigma[1L]
    sigma <- gamma_sigma[2L]
    if(sigma <= 0) return(10^6)
    if(gamma < 0 && any(X[X >= mu] >= mu - sigma/gamma)){
      return(10^6)
    }
    - sum(dgpareto(
      x = X[X >= mu],
      mu = mu,
      gamma = gamma,
      sigma = sigma,
      log = TRUE
    )) 
  }
  # estimates by method of moments
  exces <- X[X >= mu] - mu
  m <- mean(exces)
  v <- var(exces)
  scale <- m/2 * (m^2/v + 1)
  shape <- -(m^2/v - 1) / 2
  #
  inits <- c(shape, scale)
  opt <- optim(
    par = inits,
    fn = fn,
    method = "Nelder-Mead", 
    control = list(maxit = 10000L)
  )
  if(opt[["convergence"]] == 0L){
    opt[["par"]]
  }else{
    inits
  }
}

# # @importFrom stats optim
# # @noRd
# selectThreshold <- function(X, candidates){
#   params <- matrix(NA_real_, nrow = length(candidates), ncol = 2L)
#   values <- rep(NA_real_, length(candidates))
#   for(i in seq_along(candidates)){
#     mu <- candidates[i]
#     n <- length(X[X < mu])
#     p <- ifelse(n == 0L, 0, sum(X < mu)/length(X))
#     print(p)
#     pen <- ifelse(n == 0L, 0, -(length(X) -n)*(1-p) + p*n*log(mu))
#     fn <- function(gamma_sigma){
#       - sum(dgpareto(
#         x = X[X >= mu], 
#         mu = mu, 
#         gamma = gamma_sigma[1L]/(1-gamma_sigma[1L]), 
#         sigma = gamma_sigma[2L]/(1-gamma_sigma[2L]), 
#         log = TRUE
#       )) + p * log(mu) * length(X)
#     }
#     opt <- optim(
#       par = c(0.5, 0.5), fn = fn, 
#       method = "L-BFGS-B", lower = 0.01, upper = 0.99, 
#       control = list(maxit = 500L)
#     )
#     if(opt[["convergence"]] == 0L){
#       params[i, ] <- opt[["par"]]/(1-opt[["par"]])
#       values[i] <- opt[["value"]]
#     }
#   }
#   imin <- which.min(values)
#   print(values)
#   c(mu = candidates[imin], gamma = params[imin, 1L], sigma = params[imin, 2L])
# }
# 
# selectThreshold2 <- function(X, gamma, candidates){
#   sigmas <- rep(NA_real_, length(candidates))
#   values <- rep(NA_real_, length(candidates))
#   for(i in seq_along(candidates)){
#     mu <- candidates[i]
#     n <- length(X[X < mu])
#     fn <- function(sigma){
#       - sum(dgpareto(
#         x = X[X >= mu], 
#         mu = mu, 
#         gamma = gamma, 
#         sigma = sigma/(1-sigma), 
#         log = TRUE
#       )) + n*log(mu) 
#     }
#     opt <- optim(
#       par = c(0.5), fn = fn, 
#       method = "L-BFGS-B", lower = 0.01, upper = 0.99, 
#       control = list(maxit = 500L)
#     )
#     if(opt[["convergence"]] == 0L){
#       sigmas[i] <- opt[["par"]]/(1-opt[["par"]])
#       values[i] <- opt[["value"]]
#     }
#   }
#   imin <- which.min(values)
#   c(mu = candidates[imin], sigma = sigmas[imin])
# }