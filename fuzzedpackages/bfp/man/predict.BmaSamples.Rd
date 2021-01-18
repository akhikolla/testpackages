\name{predict.BmaSamples}
\alias{predict.BmaSamples}
\alias{print.predict.BmaSamples}

\title{Predict method to extract point and interval predictions from
  BmaSamples objects} 
\description{
  Predict new responses from a Bayesian model average over FP models,
  from which predictive samples have already been produced.
}
\usage{
\method{predict}{BmaSamples}(object, level=0.95, hpd=TRUE, ...)
\method{print}{predict.BmaSamples}(x, ...)
}

\arguments{
  \item{object}{valid \code{\link{BmaSamples}} object}
  \item{level}{credible level for the credible intervals (default:
    95\%)}
  \item{hpd}{should emprical hpd intervals be used (default) or simple
    quantile-based?}   
  \item{\ldots}{unused}
  \item{x}{object of S3 class \code{\link{predict.BmaSamples}}}
}

\details{
  This function summarizes the predictive samples saved in the
  \code{\link{BmaSamples}} object. Using these functions, one can obtain
  predictive credible intervals, as opposed to just using the function
  \code{\link{bmaPredict}}, which only gives the means of the predictive
  distributions.
}

\value{
  A list of class \code{predict.BmaSamples}, which has then a separate
  print method. The elements of the list are:
  \item{intervalType}{which credible intervals have been computed
    (either \dQuote{HPD} or \dQuote{equitailed})}
  \item{level}{the credible level}
  \item{newdata}{the covariate data for the predicted data points (just
    copied from \code{object})}
  \item{sampleSize}{the sample size (just copied from \code{object})}
  \item{nModels}{the number of models (just copied from \code{object})}
  \item{summaryMat}{the summary matrix for the predictions, with median,
  mean, lower and upper credible interval borders.}
}

\author{Daniel Saban\'es Bov\'e}
\seealso{\code{\link{bmaPredict}}}
\examples{
## generate a BmaSamples object
set.seed(19)

x1 <- rnorm(n=15)
x2 <- rbinom(n=15, size=20, prob=0.5) 
x3 <- rexp(n=15)

y <- rt(n=15, df=2)

test <- BayesMfp(y ~ bfp (x2, max = 4) + uc (x1 + x3), nModels = 100,
                 method="exhaustive")

## predict new responses at (again random) covariates with BMA:
testBma <- BmaSamples(test,
                      newdata=data.frame(x1 = rnorm (15),
                                         x2 = rbinom (n=15, size=5, prob=0.2) + 1,
                                         x3 = rexp (15)))
predict(testBma)
}

\keyword{methods}

