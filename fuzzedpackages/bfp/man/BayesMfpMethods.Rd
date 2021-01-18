\name{BayesMfp Methods}
\alias{BayesMfp Methods}
\alias{fitted.BayesMfp}
\alias{print.BayesMfp}
\alias{residuals.BayesMfp}

\title{Other methods for BayesMfp objects}
\description{
  Print the object (\code{print}),
  get fitted values (\code{fitted}) and corresponding residuals (\code{residuals}).
}
\usage{
\method{print}{BayesMfp}(x, ...)
\method{fitted}{BayesMfp}(object, design = getDesignMatrix(object), post =
getPosteriorParms(object, design = design), ...) 
\method{residuals}{BayesMfp}(object, ...)
}

\arguments{
  \item{x}{valid \code{\link{BayesMfp}} object}
  \item{object}{valid \code{\link{BayesMfp}} object, only the first model will be
  used.} 
  \item{design}{design matrix of the first model in the object, which
  can be supplied by the caller if it is computed beforehand}
  \item{post}{posterior parameters of the normal-gamma distribution
  (defaults to the posterior expected mean, marginalized over the
  covariance factor g)}
  \item{\ldots}{unused}
}

\author{Daniel Saban\'es Bov\'e}
\seealso{\code{\link{BayesMfp}}, \code{\link{BmaSamples Methods}}}
\examples{
## generate a BayesMfp object
set.seed(19)

x1 <- rnorm(n=15)
x2 <- rbinom(n=15, size=20, prob=0.5) 
x3 <- rexp(n=15)

y <- rt(n=15, df=2)

test <- BayesMfp(y ~ bfp (x2, max = 4) + uc (x1 + x3), nModels = 100,
                 method="exhaustive")

## the print method
test

## extract fitted values and corresponding residuals
fitted(test)
residuals(test)

}

\keyword{methods}

