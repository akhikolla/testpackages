\name{Summary of BayesMfp object}
\alias{print.summary.BayesMfp}
\alias{summary.BayesMfp}

\title{Calculate and print the summary of a BayesMfp object}
\description{
  Calculate and print the summary of a \code{\link{BayesMfp}} object,
  using S3 methods for the class. 
}
\usage{
\method{summary}{BayesMfp}(object, level=0.95, table=TRUE,
                           shrinkage=NULL, ...)
\method{print}{summary.BayesMfp}(x, ...)
}

\arguments{
  \item{object}{a valid \code{\link{BayesMfp}} object}
  \item{x}{a return value of \code{\link{summary.BayesMfp}}}
  \item{level}{credible level for coefficients HPD intervals (default: 0.95)}
  \item{table}{should a data.frame of the models be included? (default)}
  \item{shrinkage}{shrinkage factor used, where \code{NULL} defaults to
    the posterior expected shrinkage factor} 
  \item{\dots}{only used by \code{\link{summary.BayesMfp}} to pass
    arguments to \code{\link{as.data.frame.BayesMfp}}}
}
\note{
  Note that if you extract the summary of a single model with these
  functions, you ignore the uncertainty about the shrinkage factor
  t=g/(g+1) by plugging in the number \code{shrinkage}. If you want to
  incorporate this uncertainty, you must run \code{\link{BmaSamples}} on
  this model and call the corresponding method
  \code{\link{summary.BmaSamples}}. 
}
\value{
  \code{\link{summary.BayesMfp}} returns a list with S3 class
  \code{summary.BayesMfp}, where the arguments \dQuote{call},
  \dQuote{numVisited}, \dQuote{termNames}, 
  \dQuote{shiftScaleMax}, \dQuote{inclusionProbs}, \dQuote{chainlength}
  (only for model sampling results) are copied from the attributes of
  the \code{\link{BayesMfp}} object, please see its help page for
  details.  

  The other elements are:
  \item{dataframe}{the model overview as data.frame (only if
    \code{table=TRUE} was specified)}
  \item{localInclusionProbs}{local variable inclusion probability
    estimates}
  \item{nModels}{number of models contained in \code{object}}

  If there are multiple models in \code{object}, the list element
  \code{postProbs} contains the exact (for exhaustively explored model
  spaces) or estimated (if model sampling has been done) posterior model
  probabilities.

  If \code{object} contains only one FP model, then this one is
  summarized in more detail:
  \item{level}{used credible level for coefficients HPD intervals}
  \item{shrinkage}{used shrinkage factor}
  \item{summaryMat}{matrix with posterior summaries of the single
    coefficients: \dQuote{mode} gives the posterior mode,
    \dQuote{HPDlower} and \dQuote{HPDupper} give the boundaries of the HPD
    intervals with specified credible \code{level}}
  \item{sigma2Sum}{posterior summary for the regression variance: again
    mode, and lower and upper HPD bounds are given in a rowvector.}
}
\author{Daniel Saban\'es Bov\'e}
\seealso{\code{\link{summary.BmaSamples}}}
\examples{
## generate a BayesMfp object
set.seed(19)

x1 <- rnorm(n=15)
x2 <- rbinom(n=15, size=20, prob=0.5) 
x3 <- rexp(n=15)

y <- rt(n=15, df=2)

test <- BayesMfp(y ~ bfp (x2, max = 4) + uc (x1 + x3), nModels = 100,
                 method="exhaustive")

## summary of multiple models:
summary(test)

## summary of just one model (no. 10):
summary(test[10])

## internal structure is usually not interesting:
str(summary(test[10]))
}

\keyword{methods}

