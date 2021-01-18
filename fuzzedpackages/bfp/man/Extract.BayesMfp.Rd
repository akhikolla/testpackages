\name{Extract.BayesMfp}
\alias{Extract.BayesMfp}
\alias{[.BayesMfp}

\title{Extract method for BayesMfp objects}
\description{
  Extract a subset of models from a \code{\link{BayesMfp}} object.
}

\usage{
\method{[}{BayesMfp}(x, ...)
}

\arguments{
  \item{x}{valid \code{\link{BayesMfp}} object}
  \item{\dots}{transports the indexes of the models}
}

\author{Daniel Saban\'es Bov\'e}
\seealso{\code{\link{BayesMfp}}}
\examples{
## generate a BayesMfp object
set.seed(19)

x1 <- rnorm(n=15)
x2 <- rbinom(n=15, size=20, prob=0.5) 
x3 <- rexp(n=15)

y <- rt(n=15, df=2)

test <- BayesMfp(y ~ bfp (x2, max = 4) + uc (x1 + x3), nModels = 100,
                 method="exhaustive")

## extract the top ten models
test[1:10]
}

\keyword{methods}

