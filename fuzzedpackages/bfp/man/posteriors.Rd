\name{posteriors}
\alias{posteriors}

\title{Extract posterior model probability estimates from BayesMfp objects}
\description{
  Extract posterior model probability estimates (either normalized
  estimates or sampling frequencies) from \code{\link{BayesMfp}} objects.
}
\usage{
posteriors(BayesMfpObject, ind = 1)
}

\arguments{
  \item{BayesMfpObject}{a valid \code{\link{BayesMfp}} object,
    containing the models the probabilites of which one wants to
    estimate}
  \item{ind}{\code{ind = 1} means normalized posteriors, \code{ind = 2}
    means sampling frequencies}
}
\value{
  The vector of probability estimates.
}
\author{Daniel Saban\'es Bov\'e}
\examples{
## construct a BayesMfp object
set.seed(19)
x1 <- rnorm (n=15)
x2 <- rbinom (n=15, size=20, prob=0.5) 
x3 <- rexp (n=15)
y <- rt (n=15, df=2)

test <- BayesMfp (y ~ bfp (x1, max = 2) + bfp (x2, max = 2) + uc (x3), nModels = 100,
		  method="exhaustive")

## this works:
posteriors(test)

## this must not work:
## SoDA::muststop(posteriors(test, ind=2))

## only if we do model sampling there are model frequencies:
test2 <- BayesMfp (y ~ bfp (x1, max = 2) + bfp (x2, max = 2) + uc (x3), nModels = 100,
         method="sampling")
posteriors(test2, ind=2)
}
\keyword{utilities}


