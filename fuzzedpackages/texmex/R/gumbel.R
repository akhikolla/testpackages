#' @export gumbel
gumbel <- texmexFamily(name = 'Gumbel',
          log.lik = function(data, ...) {
            y <- data$y
            X.mu <- data$D$mu
            X.phi <- data$D$phi #  phi = log(sigma)

            n.mu <- ncol(X.mu)
            n.end <- n.mu + ncol(X.phi)

            function(param) {
              stopifnot(length(param) == n.end)
              mu <- X.mu %*% param[1:n.mu]
              phi <- X.phi %*% param[(1 + n.mu):n.end]
              z <- (y-mu)/exp(phi)
              sum(-exp(-z)-z-phi)
            }
          }, # Close log.lik
          param = c(mu=0, phi=0),
          info = NULL,
          sandwich = NULL,
          start = function(data){
            y <- data$y
            X.mu <- data$D$mu
            X.phi <- data$D$phi
            c(mean(y), rep(0, ncol(X.mu) - 1),
              log(IQR(y)/2), rep(0.001, ncol(X.phi)-1))
          }, # Close start

          resid = function(o){
            p <- texmexMakeParams(coef(o), o$data$D)
            (o$data$y - p[,1]) / exp(p[,2])  # Standard gumbel
          }, # Close resid

          endpoint = function(param, model){
            Inf
          },
          rl = function(m, param, model){
            param[,1] - exp(param[,2])* log(-log(1-1/m))
          }, # close rl
          delta = function(param,m,model){ # follows argument in Coles eqn (4.15) for GPD
            out <- rep(1,2) # can have vector output as call only ever assumes m length 1

            out[1] <- 1
            out[2] <- -exp(param[2]) * log(-log(1-1/m))
            out
          }, # Close delta
          density = function(n, param, model, log.d=FALSE){
            mu <- param[, 1]
            sigma <- exp(param[,2])
            dgumbel(n, mu, sigma, log.d=log.d)
          },
          rng = function(n, param, model){
            rgumbel(n, param[, 1], exp(param[, 2]))
          }, # end rng
          prob = function(n, param, model){
            pgumbel(n, param[, 1], exp(param[, 2]))
          },
          quant = function(n, param, model){
            qgumbel(n, param[, 1], exp(param[, 2]))
          }
)

#' The Gumbel distribution
#'
#' @description Density, distribution and quantile functions, and random number
#'   generation for the Gumbel distribution
#'
#' @aliases rgumbel pgumbel qgumbel dgumbel
#' @family rgumbel pgumbel qgumbel dgumbel
#'
#' @param x,q,p Vectors of quantiles or probabilities.
#' @param n The number of observations.
#' @param mu The location parameter.
#' @param sigma The scale parameter.
#' @param log.d,log.p Whether to return logged values, or to treat probabilities/densities as being logged.
#' @param lower.tail Whether to return the lower tail. If \code{lower.tail=FALSE},
#'     the upper tail is returned.
#'
#' @export
#' @name dgumbel
dgumbel <- function(x, mu, sigma, log.d=FALSE){
  xs <- -(x - mu) / sigma
  d <-  xs - exp(xs) - log(sigma)
  if (log.d){
    d
  } else {
    exp(d)
  }
}

#' @export
#' @rdname dgumbel
rgumbel <- function(n, mu, sigma){
  u <- runif(n)
  mu - sigma * log(-log(u))
}

#' @export
#' @rdname dgumbel
pgumbel <- function(q, mu, sigma, lower.tail=TRUE, log.p=FALSE){
  qs <- -(q - mu) / sigma
  if (log.p){
    res <- -exp(qs)
  } else {
    res <- exp(-exp(qs))
  }

  if (lower.tail){
    res
  } else {
    1 - res
  }
}

#' @export
#' @rdname dgumbel
qgumbel <- function(p, mu, sigma, lower.tail=TRUE, log.p = FALSE){
  if (log.p){
    p <- exp(p)
  }

  mu - sigma * log(-log(p))
}
