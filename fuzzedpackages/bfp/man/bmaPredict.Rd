\name{bmaPredict}
\alias{bmaPredict}

\title{BMA prediction for new data points}
\description{
  Make a Bayesian model averaged prediction for new data points, from
  those models saved in a \code{\link{BayesMfp}} object.
}
\usage{
bmaPredict(BayesMfpObject, postProbs = posteriors(BayesMfpObject), newdata)
}

\arguments{
  \item{BayesMfpObject}{\code{\link{BayesMfp}} object with the models
  over which the predictions should be averaged}
  \item{postProbs}{vector of posterior probabilities, which are then
  normalized to the weights of the model average (defaults to the
  normalized posterior probability estimates)}
  \item{newdata}{new covariate data as data.frame}
}
\note{
  Note that this function is not an S3 predict method for
  \code{\link{BmaSamples}} objects, but a function working on
  \code{\link{BayesMfp}} objects (because we do not need BMA samples to
  do BMA point predictions).
}
\value{
  The predicted values as a vector.
}
\author{Daniel Saban\'es Bov\'e}
\seealso{\code{\link{BmaSamples Methods}}}
\examples{
## generate a BayesMfp object
set.seed(19)

x1 <- rnorm(n=15)
x2 <- rbinom(n=15, size=20, prob=0.5) 
x3 <- rexp(n=15)

y <- rt(n=15, df=2)

test <- BayesMfp(y ~ bfp (x2, max = 4) + uc (x1 + x3), nModels = 100,
                 method="exhaustive")

## predict new responses at (again random) covariates
bmaPredict(test,
           newdata = list(x1 = rnorm(n=15),
                          x2 = rbinom(n=15, size=5, prob=0.2) + 1,
                          x3 = rexp(n=15)))
}
\keyword{regression}

