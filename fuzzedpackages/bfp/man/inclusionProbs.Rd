\name{inclusionProbs}
\alias{inclusionProbs}

\title{ Compute (model averaged) posterior variable inclusion probabilites}
\description{
  Compute (model averaged) posterior inclusion probabilites for the
  uncertain variables (including FP variables) based on a
  \code{\link{BayesMfp}} object. 
}
\usage{
inclusionProbs(BayesMfpObject, postProbs = posteriors(BayesMfpObject, ind = 1))
}

\arguments{
  \item{BayesMfpObject}{valid \code{\link{BayesMfp}} object}
  \item{postProbs}{posterior probabilities to weight the models
    (defaults to the normalized probability estimates)}
}
\value{
  Named numeric vector with the estimated variable inclusion
  probabilities. Note that these can differ noticeably from the
  \dQuote{global} inclusion probabilities computed from all discovered
  models, from which only the best were saved in the
  \code{\link{BayesMfp}} object.
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

## now get the local inclusion probabilities
local <- inclusionProbs(test)

## they can be compared with the global inclusion probabilities
local - attr(test, "inclusionProbs")
}

\keyword{htest}
