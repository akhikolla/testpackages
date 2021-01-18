\name{dinvGamma}
\alias{dinvGamma}
\alias{pinvGamma}
\alias{qinvGamma}
\alias{rinvGamma}

\title{Functions for the inverse gamma distribution.}
\description{
  Functions for the inverse gamma distribution, analogues to the gamma
  distribution functions.
}
\usage{
dinvGamma(x, a, b, log = FALSE, normalize = TRUE)
pinvGamma(q, a, b, lower.tail = TRUE, log.p = FALSE)
qinvGamma(p, a, b, lower.tail = TRUE, log.p = FALSE)
rinvGamma(n, a, b)
}

\arguments{
  \item{x, q}{vector of quantiles}
  \item{a}{shape parameter}
  \item{b}{rate parameter}
  \item{log, log.p}{logical; if \code{TRUE}, probabilities/densities p
    are returned as log(p)} 
  \item{normalize}{normalize the density function? (default)}
  \item{lower.tail}{logical; if \code{TRUE} (default), probabilities are P[X <= x],
          otherwise, P[X > x]}
  \item{p}{vector of probabilities}
  \item{n}{number of observations. If \code{length(n) > 1}, the length is
          taken to be the number required}
}
\details{
  See the gamma distribution (\code{\link{rgamma}}) for the details.
}
\author{Daniel Saban\'es Bov\'e}
\keyword{distribution}
\keyword{internal}
