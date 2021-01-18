warpSpline2 <- list(
  fun = function(z, par, L = 2) {
    if (length(par) != L) stop("'par' must be of length ", L)
    xknots <- seq(from = 0, to = 1, length.out = L)
    if (!requireNamespace("DiceKriging", quietly = TRUE)) {
      stop("DiceKriging required")
    }
    y <- DiceKriging::scalingFun1d(z, knots = xknots, eta = par)
    attr(y, "gradient") <- DiceKriging:::scalingGrad1d(z, knots = xknots)
    return(y)
  },
  parNames = paste("eta"),
  parDefault = c("eta" = 1),
  parLower = c("eta" = 1e-10),
  parUpper = c("eta" = Inf),
  hasGrad = TRUE,
  isCdf = FALSE
)

warpSpline1 <- list(
  fun = function(z, par, L = 2) {
    if (length(par) != (L-1)) stop("'par' must be of length ", L-1)
    xknots <- seq(from = 0, to = 1, length.out = L)
    yknots <- c(0, cumsum(par))
    y <- approx(x = xknots, y = yknots, xout = z)$y
    attr(y, "gradient") <- outer(z, xknots[-1], function(z, t){z >= t}) * 1
    return(y)
  },
  parNames = paste("eta"),
  parDefault = c("eta" = 1),
  parLower = c("eta" = 1e-10),
  parUpper = c("eta" = Inf),
  hasGrad = TRUE,
  isCdf = FALSE
)

warpPower <- list(
  fun = function(z, par, L) {
    y <- pbeta(q = z, shape1 = par[1], shape2 = 1)
    ind <- (z > 0) & (z < 1)
    grad <- rep(0, length(z))
    if (any(ind)) grad[ind] <- y[ind] * log(z[ind])
    attr(y, "gradient") <- matrix(grad, ncol = 1)
    return(y)
  },
  parNames = "pow",
  parDefault = c("pow" = 1),
  parLower = c("pow" = 1e-10),
  parUpper = c("pow" = Inf),
  hasGrad = TRUE,
  isCdf = TRUE
)

eps <- 1e-10
warpNorm <- list(
  fun = function(z, par, L) {
    Az <- pnorm(z, mean = par[1], sd = par[2])
    A1 <- pnorm(1, mean = par[1], sd = par[2])
    A0 <- pnorm(0, mean = par[1], sd = par[2])
    az <- dnorm(z, mean = par[1], sd = par[2])
    a1 <- dnorm(1, mean = par[1], sd = par[2])
    a0 <- dnorm(0, mean = par[1], sd = par[2])
    N <- Az - A0
    D <- A1 - A0
    y <- (Az - A0)/D
    grad <- matrix(0, nrow = length(z), ncol = 2)
    colnames(grad) <- c("mean", "sd")
    grad[, 1] <-  (az - a0) - N * (a1 - a0) / D
    grad[, 1] <- - grad[, 1] / D

    grad[, 2] <-  (az * (z - par[1]) - a0 * (0 - par[1])) - 
                  N * (a1 * (1 - par[1]) - a0 * (0 - par[1])) / D
    grad[, 2] <- - grad[, 2] / par[2] / D
    attr(y, "gradient") <- grad
    return(y)
  },
  parNames = c("mean", "sd"),
  parDefault = c(mean = 0.5, sd = 3),
  parLower = c(mean = eps, sd = eps),
  parUpper = c(mean = 1 - eps, sd = Inf),
  hasGrad = TRUE,
  isCdf = TRUE
)

warpUnorm <- list(
  fun = function(z, par, L) {
    y <- pnorm(z, mean = par[1], sd = par[2])
    grad <- matrix(0, nrow = length(z), ncol = 2)
    colnames(grad) <- c("mean", "sd")
    phi <- dnorm(z, mean = par[1], sd = par[2])
    grad[, 1] <- - phi 
    grad[, 2] <- - phi * (z - par[1]) / par[2]
    attr(y, "gradient") <- grad
    return(y)
  },
  parNames = c("mean", "sd"),
  parDefault = c(mean = 0.5, sd = 3),
  parLower = c(mean = eps, sd = eps),
  parUpper = c(mean = 1 - eps, sd = Inf),
  hasGrad = TRUE,
  isCdf = FALSE
)





#class(warpPower) <- c("warp", "list")


# warpNorm <- list(
#   fun = function(z, par) {
#     y <- pnorm(q = z, mean = par[1], sd = par[2])
#     f <- dnorm(q = z, mean = par[1], sd = par[2])
#     ind <- (z >= 0) & (z <= 1)
#     grad <- matrix(0, nrow = length(z), ncol = 2)
#     colnames(grad) <- c("mean", "sd")
#     if (any(ind)) {
#       grad[ind, 1] <- - f[ind] / par[2]
#       grad[ind, 2] <- - f[ind] * (z[ind] - par[1]) / par[2]^2
#     }
#     attr(y, "gradient") <- grad
#     return(y)
#   },
#   parNames = "pow",
#   parDefault = c(pow = 1),
#   parLower = c(pow = 1e-10),
#   parUpper = c(pow = Inf),
#   hasGrad = TRUE
# )
#class(warpPower) <- c("warp", "list")
