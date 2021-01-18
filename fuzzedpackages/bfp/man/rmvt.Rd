\name{rmvt}
\alias{rmvt}

\title{Multivariate Student Random Deviates}
\description{
  This function provides a random vector generator for the multivariate
  Student distribution with mean vector \code{mean}, scale matrix
  \code{sigma} and degrees of freedom \code{df}. 
}
\usage{
rmvt(n, sigma = diag(2), mu = rep(0, 2), df = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{Number of observations.}
  \item{sigma}{symmetric and positive definite scale matrix}
  \item{mu}{mean vector}
  \item{df}{degree of freedom as integer}
}
\value{
  The vector samples, stacked in the rows of a matrix.
}
\note{
  Note that the covariance matrix is only defined for \code{df > 2}, and
  is given by \code{sigma * df / (df - 2)}.
}

\author{Daniel Saban\'es Bov\'e}
\examples{
## samples from the multivariate Cauchy distribution:
bfp:::rmvt(20)

## here the covariance exists:
sigma <- matrix(c(1, 0.5, 0.5, 1), nrow=2)
df <- 10

## theoretical covariance:
sigma * df / (df - 2)

## this should be close:
cov(bfp:::rmvt(n=100000, sigma=sigma, df=df))
}
\keyword{distribution}
\keyword{internal}
