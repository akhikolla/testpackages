\name{getPosteriorParms}
\alias{getPosteriorParms}

\title{Extract updated posterior parameters for the normal inverse gamma
  distribution from a model, given a shrinkage factor.} 
\description{
  Conditional on a fixed shrinkage factor t=g/(g+1), the posterior joint distribution
  of the effects and the regression variance is normal inverse gamma.
  With this function, you can compute the parameters of this distribution.
}
\usage{
getPosteriorParms(x, shrinkage=x[[1]]$postExpectedShrinkage,  
                  design = getDesignMatrix(x))
}

\arguments{
  \item{x}{a valid \code{\link{BayesMfp}}-Object, only first list element will be recognized}
  \item{shrinkage}{shrinkage factor used in the computations (defaults
  to the posterior expected shrinkage factor in the model \code{x[1]})} 
  \item{design}{(centered) design matrix for the model}
}
\value{
  A list with four parameters:
  \item{aStar}{the first parameter of the inverse gamma distribution}
  \item{VStar}{the covariance matrix part of the multivariate normal
    distribution}
  \item{mStar}{the expectation of the multivariate normal
    distribution}
  \item{bStar}{the second parameter of the inverse gamma distribution}
}
\author{Daniel Saban\'es Bov\'e}
\examples{
## construct a BayesMfp object
set.seed(19)

x1 <- rnorm (n=15)
x2 <- rbinom (n=15, size=20, prob=0.5) 
x3 <- rexp (n=15)

y <- rt (n=15, df=2)

test <- BayesMfp (y ~ bfp (x2, max = 4) + uc (x1 + x3), nModels = 200, method="exhaustive")

## now get the posterior parameters of the third best model
getPosteriorParms(test[3])
}
\keyword{regression}

