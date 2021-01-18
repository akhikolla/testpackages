\name{Summary of BmaSamples object}
\alias{print.summary.BmaSamples}
\alias{summary.BmaSamples}

\title{Calculate and print the summary of a BmaSamples object}
\description{
  Calculate and print the summary of a \code{\link{BmaSamples}} object,
  using S3 methods for the class. 
}
\usage{
\method{summary}{BmaSamples}(object, level = 0.95, hpd = TRUE, ...)
\method{print}{summary.BmaSamples}(x, table = TRUE, ...)
}

\arguments{
  \item{object}{a valid \code{\link{BmaSamples}} object}
  \item{level}{credible level for coefficients credible intervals}
  \item{hpd}{should emprical hpd intervals be used (default) or simple quantile-based?}
  \item{x}{a return value of \code{\link{summary.BmaSamples}}}
  \item{table}{should the model table been shown? (default)}
  \item{\dots}{unused}
}
\value{
  The summary method returns an S3 object, where \dQuote{sampleSize},
  \dQuote{modelData} and \dQuote{modelFreqs} are copied from the
  \code{\link{BmaSamples}} object, please see its help page for the
  details. \dQuote{intervalType} and \dQuote{level} copy the function's
  parameters. 

  \dQuote{summaryMat} contains the posterior summaries for the intercept
  and uncertain fixed form covariates. \dQuote{sigma2Sum} and
  \dQuote{shrinkageSum} contain the posterior summaries for the
  regression variance and the shrinkage factor, respectively. The
  summaries are always the median, mean, lower and upper credible
  bounds for the coefficients.
}
\author{Daniel Saban\'es Bov\'e}
\seealso{\code{\link{summary.BayesMfp}}}
\examples{
## generate a BmaSamples object
set.seed(19)

x1 <- rnorm(n=15)
x2 <- rbinom(n=15, size=20, prob=0.5) 
x3 <- rexp(n=15)

y <- rt(n=15, df=2)

test <- BayesMfp(y ~ bfp (x2, max = 4) + uc (x1 + x3), nModels = 100,
                 method="exhaustive")

testBma <- BmaSamples(test)

## look at the summary
summary(testBma)

## and its structure
str(summary(testBma))
}

\keyword{methods}

