# ================================ GEV functions ===============================

#' The Generalised Extreme Value Distribution
#'
#' Density function, distribution function, quantile function and
#' random generation for the generalised extreme value (GEV)
#' distribution.
#'
#' @param x,q Numeric vectors of quantiles.
#' @param p A numeric vector of probabilities in [0,1].
#' @param loc,scale,shape Numeric vectors.
#'   Location, scale and shape parameters.
#'   All elements of \code{scale} must be positive.
#' @param n Numeric scalar.  The number of observations to be simulated.
#'   If \code{length(n) > 1} then \code{length(n)} is taken to be the number
#'   required.
#' @param log,log.p A logical scalar; if TRUE, probabilities p are given as
#'   log(p).
#' @param lower.tail A logical scalar.  If TRUE (default), probabilities
#'   are \eqn{P[X \leq x]}{P[X <= x]}, otherwise, \eqn{P[X > x]}{P[X > x]}.
#' @param m A numeric scalar.  The distribution is reparameterised by working
#'  with the GEV(\code{loc, scale, shape}) distribution function raised to the
#'  power \code{m}.  See \strong{Details}.
#' @details The distribution function of a GEV distribution with parameters
#'  \code{loc} = \eqn{\mu}, \code{scale} = \eqn{\sigma (> 0)} and
#'  \code{shape} = \eqn{\xi} is
#'  \deqn{F(x) = \exp\{-[1 + \xi (x - \mu) / \sigma] ^ {-1/\xi} \}}{%
#'        F(x) = exp{ -[1 + \xi (x - \mu) / \sigma] ^ (-1/\xi)} }
#'  for \eqn{1 + \xi (x - \mu) / \sigma > 0}.  If \eqn{\xi = 0} the
#'  distribution function is defined as the limit as \eqn{\xi} tends to zero.
#'  The support of the distribution depends on \eqn{\xi}: it is
#'  \eqn{x \leq \mu - \sigma / \xi}{x <= \mu - \sigma / \xi} for \eqn{\xi < 0};
#'  \eqn{x \geq \mu - \sigma / \xi}{x >= \mu - \sigma / \xi} for \eqn{\xi > 0};
#'  and \eqn{x} is unbounded for \eqn{\xi = 0}.
#'  Note that if \eqn{\xi < -1} the GEV density function becomes infinite
#'  as \eqn{x} approaches \eqn{\mu -\sigma / \xi} from below.
#'
#'  If \code{lower.tail = TRUE} then if \code{p = 0} (\code{p = 1}) then
#'  the lower (upper) limit of the distribution is returned, which is
#'  \code{-Inf} or \code{Inf} in some cases.  Similarly, but reversed,
#'  if \code{lower.tail = FALSE}.
#'
#'  See
#'  \url{https://en.wikipedia.org/wiki/Generalized_extreme_value_distribution}
#'  for further information.
#'
#'  The effect of \code{m} is to change the location, scale and shape
#'  parameters to
#'  \eqn{(\mu + \sigma \log m, \sigma, \xi)}{(\mu + \sigma log m, \sigma, \xi)}
#'  if \eqn{\xi = 0} and
#'  \eqn{(\mu + \sigma (m ^ \xi - 1) / \xi, \sigma m ^ \xi, \xi)}.
#'  For integer \code{m} we can think of this as working with the
#'  maximum of \code{m} independent copies of the original
#'  GEV(\code{loc, scale, shape}) variable.
#' @return \code{dgev} gives the density function, \code{pgev} gives the
#'   distribution function, \code{qgev} gives the quantile function,
#'   and \code{rgev} generates random deviates.
#'
#'   The length of the result is determined by \code{n} for \code{rgev},
#'   and is the maximum of the lengths of the numerical arguments for the
#'   other functions.
#'
#'   The numerical arguments other than \code{n} are recycled to the length
#'   of the result.
#' @references Jenkinson, A. F. (1955) The frequency distribution of the
#'   annual maximum (or minimum) of meteorological elements.
#'   \emph{Quart. J. R. Met. Soc.}, \strong{81}, 158-171.
#'   Chapter 3: \url{https://doi.org/10.1002/qj.49708134804}
#' @references Coles, S. G. (2001) \emph{An Introduction to Statistical
#'   Modeling of Extreme Values}, Springer-Verlag, London.
#'   \url{https://doi.org/10.1007/978-1-4471-3675-0_3}
#' @examples
#' dgev(-1:4, 1, 0.5, 0.8)
#' dgev(1:6, 1, 0.5, -0.2, log = TRUE)
#' dgev(1, shape = c(-0.2, 0.4))
#'
#' pgev(-1:4, 1, 0.5, 0.8)
#' pgev(1:6, 1, 0.5, -0.2)
#' pgev(1, c(1, 2), c(1, 2), c(-0.2, 0.4))
#' pgev(-3, c(1, 2), c(1, 2), c(-0.2, 0.4))
#' pgev(7, 1, 1, c(-0.2, 0.4))
#'
#' qgev((1:9)/10, 2, 0.5, 0.8)
#' qgev(0.5, c(1,2), c(0.5, 1), c(-0.5, 0.5))
#'
#' p <- (1:9)/10
#' pgev(qgev(p, 1, 2, 0.8), 1, 2, 0.8)
#'
#' rgev(6, 1, 0.5, 0.8)
#' @name gev
NULL
## NULL

# ----------------------------- dgev ---------------------------------

#' @rdname gev
#' @export
dgev <- function (x, loc = 0, scale = 1, shape = 0, log = FALSE, m = 1) {
  if (any(scale <= 0)) {
    stop("invalid scale: scale must be positive.")
  }
  if (length(x) == 0) {
    return(numeric(0))
  }
  max_len <- max(length(x), length(loc), length(scale), length(shape),
                 length(m))
  x <- rep_len(x, max_len)
  loc <- rep_len(loc, max_len)
  scale <- rep_len(scale, max_len)
  shape <- rep_len(shape, max_len)
  m <- rep_len(m, max_len)
  x <- (x - loc) / scale
  xx <- 1 + shape * x
  d <- ifelse(xx < 0 | is.infinite(x), -Inf,
         ifelse(xx == 0 & shape == -1, 0,
           ifelse(xx == 0 & shape < -1, Inf,
             ifelse(xx == 0, -Inf,
               ifelse(abs(shape) > 1e-6,
                 -(1 + 1 / shape) * logNegNA(xx) - m * xx ^ (-1/ shape),
      -x + shape * x * (x - 2) / 2 - m * exp(-x + shape * x ^ 2 / 2))))))
  d <- d + log(m) - log(scale)
  if (!log) {
    d <- exp(d)
  }
  return(d)
}

# ----------------------------- pgev ---------------------------------

#' @rdname gev
#' @export
pgev <- function(q, loc = 0, scale = 1, shape = 0, lower.tail = TRUE,
                 log.p = FALSE, m = 1) {
  if (any(scale <= 0)) {
    stop("invalid scale: scale must be positive.")
  }
  if (length(q) == 0) {
    return(numeric(0))
  }
  max_len <- max(length(q), length(loc), length(scale), length(shape),
                 length(m))
  q <- rep_len(q, max_len)
  loc <- rep_len(loc, max_len)
  scale <- rep_len(scale, max_len)
  shape <- rep_len(shape, max_len)
  m <- rep_len(m, max_len)
  q <- (q - loc) / scale
  p <- 1 + shape * q
  p <- ifelse(abs(shape) > 1e-6,
              -m * pmax(p, 0) ^ (-1 / shape),
              ifelse(is.infinite(q), log((1 + sign(q)) / 2),
                     -m * exp(-q + shape * q ^ 2 / 2)))
  if (lower.tail) {
    if (!log.p) {
      p <- exp(p)
    }
  } else {
    p <- -expm1(p)
    if (log.p) {
      p <- log(p)
    }
  }
  return(p)
}

# ----------------------------- qgev ---------------------------------

#' @rdname gev
#' @export
qgev <- function (p, loc = 0, scale = 1, shape = 0, lower.tail = TRUE,
                  log.p = FALSE, m = 1) {
  if (any(scale <= 0)) {
    stop("invalid scale: scale must be positive.")
  }
  if (length(p) == 0) {
    return(numeric(0))
  }
  if (!log.p & any(p < 0 | p > 1, na.rm = TRUE)) {
    stop("invalid p: p must be in [0,1].")
  }
  max_len <- max(length(p), length(loc), length(scale), length(shape),
                 length(m))
  p <- rep_len(p, max_len)
  loc <- rep_len(loc, max_len)
  scale <- rep_len(scale, max_len)
  shape <- rep_len(shape, max_len)
  m <- rep_len(m, max_len)
  if (log.p) {
    p <- exp(p)
  }
  if (!lower.tail) {
    p <- 1 - p
  }
  xp <- -log(p) / m
  mult <- box_cox_vec(x = xp, lambda = -shape)
  return(loc - scale * mult)
}

# ----------------------------- rgev ---------------------------------

#' @rdname gev
#' @export
rgev <- function (n, loc = 0, scale = 1, shape = 0, m = 1) {
  max_len <- ifelse(length(n) > 1, length(n), n)
  loc <- rep_len(loc, max_len)
  scale <- rep_len(scale, max_len)
  shape <- rep_len(shape, max_len)
  m <- rep_len(m, max_len)
  return(qgev(stats::runif(n), loc = loc, scale = scale, shape = shape, m = m))
}

# ================================ GP functions ===============================

#' The Generalised Pareto Distribution
#'
#' Density function, distribution function, quantile function and
#' random generation for the generalised Pareto (GP) distribution.
#'
#' @param x,q Numeric vectors of quantiles.  All elements of \code{x}
#'   and \code{q} must be non-negative.
#' @param p A numeric vector of probabilities in [0,1].
#' @param loc,scale,shape Numeric vectors.
#'   Location, scale and shape parameters.
#'   All elements of \code{scale} must be positive.
#' @param n Numeric scalar.  The number of observations to be simulated.
#'   If \code{length(n) > 1} then \code{length(n)} is taken to be the number
#'   required.
#' @param log,log.p A logical scalar; if TRUE, probabilities p are given as
#'   log(p).
#' @param lower.tail A logical scalar.  If TRUE (default), probabilities
#'   are \eqn{P[X \leq x]}{P[X <= x]}, otherwise, \eqn{P[X > x]}{P[X > x]}.
#' @details The distribution function of a GP distribution with parameters
#'  \code{location} = \eqn{\mu}, \code{scale} = \eqn{\sigma (> 0)} and
#'  \code{shape} = \eqn{\xi} is
#'  \deqn{F(x) = 1 - [1 + \xi (x - \mu) / \sigma] ^ {-1/\xi}}{%
#'        F(x) = 1 - [1 + \xi (x - \mu) / \sigma] ^ (-1/\xi)}
#'  for \eqn{1 + \xi (x - \mu) / \sigma > 0}.  If \eqn{\xi = 0} the
#'  distribution function is defined as the limit as \eqn{\xi} tends to zero.
#'  The support of the distribution depends on \eqn{\xi}: it is
#'  \eqn{x \geq \mu}{x >= \mu} for \eqn{\xi \geq 0}{\xi >= 0}; and
#'  \eqn{\mu \leq x \leq \mu - \sigma / \xi}{\mu <= x <= \mu - \sigma / \xi}
#'  for \eqn{\xi < 0}.  Note that if \eqn{\xi < -1} the GP density function
#'  becomes infinite as \eqn{x} approaches \eqn{\mu - \sigma/\xi}.
#'
#'  If \code{lower.tail = TRUE} then if \code{p = 0} (\code{p = 1}) then
#'  the lower (upper) limit of the distribution is returned.
#'  The upper limit is \code{Inf} if \code{shape} is non-negative.
#'  Similarly, but reversed, if \code{lower.tail = FALSE}.
#'
#'  See
#'  \url{https://en.wikipedia.org/wiki/Generalized_Pareto_distribution}
#'   for further information.
#' @return \code{dgp} gives the density function, \code{pgp} gives the
#'   distribution function, \code{qgp} gives the quantile function,
#'   and \code{rgp} generates random deviates.
#' @references Pickands, J. (1975) Statistical inference using extreme
#'   order statistics. \emph{Annals of Statistics}, \strong{3}, 119-131.
#'   \url{https://doi.org/10.1214/aos/1176343003}
#' @references Coles, S. G. (2001) \emph{An Introduction to Statistical
#'   Modeling of Extreme Values}, Springer-Verlag, London.
#'   Chapter 4: \url{https://doi.org/10.1007/978-1-4471-3675-0_4}
#' @examples
#' dgp(0:4, scale = 0.5, shape = 0.8)
#' dgp(1:6, scale = 0.5, shape = -0.2, log = TRUE)
#' dgp(1, scale = 1, shape = c(-0.2, 0.4))
#'
#' pgp(0:4, scale = 0.5, shape = 0.8)
#' pgp(1:6, scale = 0.5, shape = -0.2)
#' pgp(1, scale = c(1, 2), shape = c(-0.2, 0.4))
#' pgp(7, scale = 1, shape = c(-0.2, 0.4))
#'
#' qgp((0:9)/10, scale = 0.5, shape = 0.8)
#' qgp(0.5, scale = c(0.5, 1), shape = c(-0.5, 0.5))
#'
#' p <- (1:9)/10
#' pgp(qgp(p, scale = 2, shape = 0.8), scale = 2, shape = 0.8)
#'
#' rgp(6, scale = 0.5, shape = 0.8)
#' @name gp
NULL
## NULL

# ----------------------------- dgp ---------------------------------

#' @rdname gp
#' @export
dgp <- function (x, loc = 0, scale = 1, shape = 0, log = FALSE) {
  if (any(scale < 0)) {
    stop("invalid scale: scale must be positive.")
  }
  if (length(x) == 0) {
    return(numeric(0))
  }
  max_len <- max(length(x), length(loc), length(scale), length(shape))
  x <- rep_len(x, max_len)
  loc <- rep_len(loc, max_len)
  scale <- rep_len(scale, max_len)
  shape <- rep_len(shape, max_len)
  x <- (x - loc) / scale
  xx <- 1 + shape * x
  x <- ifelse(x < 0 | xx < 0 | is.infinite(x), -Inf,
              ifelse(xx == 0 & shape == -1, 0,
                     ifelse(xx == 0 & shape < -1, Inf,
                            ifelse(xx == 0, -Inf,
                                   ifelse(abs(shape) > 1e-6,
                                          -(1 + 1 / shape) * logNegNA(xx),
                                          -x + shape * x * (x - 2) / 2)))))
  x <- x - log(scale)
  if (!log) {
    x <- exp(x)
  }
  return(x)
}

# ----------------------------- pgp ---------------------------------

#' @rdname gp
#' @export
pgp <- function (q, loc = 0, scale = 1, shape = 0, lower.tail = TRUE,
                 log.p = FALSE){
  if (any(scale < 0)) {
    stop("invalid scale: scale must be positive.")
  }
  if (length(q) == 0) {
    return(numeric(0))
  }
  max_len <- max(length(q), length(loc), length(scale), length(shape))
  q <- rep_len(q, max_len)
  loc <- rep_len(loc, max_len)
  scale <- rep_len(scale, max_len)
  shape <- rep_len(shape, max_len)
  q <- pmax(q - loc, 0) / scale
  p <- 1 + shape * q
  p <- ifelse(abs(shape) > 1e-6,
              1 - pmax(p, 0) ^ (-1 / shape),
              ifelse(is.infinite(q), (1 + sign(q)) / 2,
                     1 - exp(-q + shape * q ^ 2 / 2)))
  if (lower.tail) {
    if (log.p) {
      p <- log(p)
    }
  } else {
    if (log.p) {
      p <- log1p(-p)
    } else {
      p <- 1 - p
    }
  }
  return(p)
}

# ----------------------------- qgp ---------------------------------

#' @rdname gp
#' @export
qgp <- function (p, loc = 0, scale = 1, shape = 0, lower.tail = TRUE,
                 log.p = FALSE) {
  if (any(scale < 0)) {
    stop("invalid scale: scale must be positive.")
  }
  if (length(p) == 0) {
    return(numeric(0))
  }
  if (!log.p & any(p < 0 | p > 1, na.rm = TRUE)) {
    stop("invalid p: p must be in [0,1].")
  }
  max_len <- max(length(p), length(loc), length(scale), length(shape))
  p <- rep_len(p, max_len)
  loc <- rep_len(loc, max_len)
  scale <- rep_len(scale, max_len)
  shape <- rep_len(shape, max_len)
  if (log.p) {
    p <- exp(p)
  }
  if (!lower.tail) {
    p <- 1 - p
  }
  mult <- box_cox_vec(x = 1 - p, lambda = -shape)
  return(loc - scale * mult)
}

# ----------------------------- rgp ---------------------------------

#' @rdname gp
#' @export
rgp <- function (n, loc = 0, scale = 1, shape = 0) {
  max_len <- ifelse(length(n) > 1, length(n), n)
  loc <- rep_len(loc, max_len)
  scale <- rep_len(scale, max_len)
  shape <- rep_len(shape, max_len)
  return(qgp(stats::runif(n), loc = loc, scale = scale, shape = shape))
}

# =========================== binomial-GP functions ============================

# --------------------------- dbingp ---------------------------------

dbingp <- function(x, p_u = 0.5 , loc = 0, scale = 1, shape = 0, log = FALSE) {
  #
  # Binomial-GP `density' function.
  #`
  # Args:
  #   x     : Numeric vector of quantiles.  No element of x can be < loc.
  #   p_u   : Numeric vector of threshold exceedance probabilities in (0,1).
  #   loc   : Numeric vector of GP location parameters: usually the threshold.
  #   scale : Numeric vector of GP scale parameters.
  #   shape : Numeric vector of GP shape parameters.
  #   log   : A logical scalar.  If TRUE the log-density is returned.
  #
  if (any(x < loc)) {
    stop("Invalid x: no element of  can be less than loc.")
  }
  if (any(scale < 0)) {
    stop("invalid scale: scale must be positive.")
  }
  if (any(p_u <= 0) || any(p_u >= 1)) {
    stop("invalid p_u: p_u must be in (0,1).")
  }
  d <- dgp(x = x, loc = loc, scale = scale, shape = shape, log = TRUE) +
    log(p_u)
  if (!log) {
    d <- exp(d)
  }
  return(d)
}

# --------------------------- pbingp ---------------------------------

pbingp <- function(q, p_u = 0.5 , loc = 0, scale = 1, shape = 0,
                   lower.tail = TRUE) {
  #
  # Binomial-GP distribution function.
  #`
  # Args:
  #   q          : Numeric vector of quantiles.  No element of q can be < loc.
  #   p_u        : Numeric vector of threshold exceedance probabilities in
  #                (0,1).
  #   loc        : Numeric vector of GP location parameters: usually the
  #               threshold.
  #   scale      : Numeric vector of GP scale parameters.
  #   shape      : Numeric vector of GP shape parameters.
  #   lower.tail : A logical scalar.  If TRUE (default), probabilities
  #                are P[X <= x], otherwise, P[X > x].
  #
  if (any(q < loc)) {
    stop("Invalid q: no element of q can be less than loc.")
  }
  if (min(scale) < 0) {
    stop("invalid scale: scale must be positive.")
  }
  if (any(p_u <= 0) || any(p_u >= 1)) {
    stop("invalid p_u: p_u must be in (0,1).")
  }
  p <- 1 - p_u * pgp(q = q, loc = loc, scale = scale, shape = shape,
                     lower.tail = FALSE)
  if (!lower.tail) {
    p <- 1 - p
  }
  return(p)
}

# --------------------------- qbingp ---------------------------------

qbingp <- function(p, p_u = 0.5, loc = 0, scale = 1, shape = 0,
                   lower.tail = TRUE) {
  #
  # Binomial-GP quantiles.
  #`
  # Args:
  #   p          : Numeric vector of probabilities in [0,1].
  #   p_u        : Numeric vector of threshold exceedance probabilities in
  #                (0,1].
  #   loc        : Numeric vector of GP location parameters: usually the
  #                threshold.
  #   scale      : Numeric vector of GP scale parameters.
  #   shape      : Numeric vector of GP shape parameters.
  #   lower.tail : A logical scalar.  If TRUE (default), probabilities
  #                are P[X <= x], otherwise, P[X > x].
  #
  if (!lower.tail) {
    p <- 1 - p
  }
  if (any(p < 1 - p_u)) {
    if (lower.tail) {
      stop("Invalid p: no element of p can be less than 1 - p_u.")
    }
    if (!lower.tail) {
      stop("Invalid p: no element of p can be greater than than p_u.")
    }
  }
  if (min(scale) < 0) {
    stop("invalid scale: scale must be positive.")
  }
  if (any(p_u <= 0) || any(p_u > 1)) {
    stop("invalid p_u: p_u must be in (0,1].")
  }
  pnew <- 1 - (1 - p) / p_u
  q <- qgp(p = pnew, loc = loc, scale = scale, shape = shape)
  return(q)
}
