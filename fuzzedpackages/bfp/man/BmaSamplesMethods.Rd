\name{BmaSamples Methods}
\alias{BmaSamples Methods}
\alias{fitted.BmaSamples}
\alias{print.BmaSamples}
\alias{residuals.BmaSamples}

\title{Other methods for BmaSamples objects}
\description{
  Print the object (\code{print}),
  get fitted values (\code{fitted}) and corresponding residuals
  (\code{residuals}).
}
\usage{
\method{print}{BmaSamples}(x, ...)
\method{fitted}{BmaSamples}(object, ...)
\method{residuals}{BmaSamples}(object, ...)
}

\arguments{
  \item{x}{valid \code{\link{BmaSamples}} object}
  \item{object}{valid \code{\link{BmaSamples}} object}
  \item{\dots}{unused}
}
\author{Daniel Saban\'es Bov\'e}
\seealso{\code{\link{predict.BmaSamples}}, \code{\link{summary.BmaSamples}}}
\examples{
## construct a BayesMfp object
set.seed(19)

x1 <- rnorm (n=15)
x2 <- rbinom (n=15, size=20, prob=0.5) 
x3 <- rexp (n=15)

y <- rt (n=15, df=2)

test <- BayesMfp (y ~ bfp (x2, max = 4) + uc (x1 + x3), nModels = 200, method="exhaustive")

## now draw samples from the Bayesian model average
testBma <- BmaSamples (test)

## the print method:
testBma

## the fitted method:
fitted(testBma)

## and the corresponding residuals:
residuals(testBma)

}

\keyword{methods}
