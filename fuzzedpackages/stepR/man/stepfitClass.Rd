\name{stepfit}
\alias{[.stepfit}
\alias{stepfit}
\alias{print.stepfit}
\alias{plot.stepfit}
\alias{lines.stepfit}
\alias{fitted.stepfit}
\alias{residuals.stepfit}
\alias{logLik.stepfit}

\title{Fitted step function}

\description{
Constructs an object containing a step function fitted to some data.
}

\usage{
stepfit(cost, family, value, param = NULL, leftEnd, rightEnd, x0,
  leftIndex = leftEnd, rightIndex = rightEnd)
\method{[}{stepfit}(x, i, j, drop = if(missing(i)) TRUE else
  if(missing(j)) FALSE else length(j) == 1, refit = FALSE) 
\method{print}{stepfit}(x, \dots)
\method{plot}{stepfit}(x, dataspace = TRUE, \dots)
\method{lines}{stepfit}(x, dataspace = TRUE, \dots)
\method{fitted}{stepfit}(object, \dots)
\method{residuals}{stepfit}(object, y, \dots)
\method{logLik}{stepfit}(object, df = NULL, nobs = object$rightIndex[nrow(object)], \dots)
}

\arguments{
\item{cost}{the value of the cost-functional used for the fit: RSS for family \code{gauss}, log-likelihood (up to a constant) for families \code{poisson} and \code{binomial}}
\item{family}{distribution of the errors, either \code{"gauss"}, \code{"poisson"} or \code{"binomial"}}
\item{value}{a numeric vector containing the fitted values for each block; its length gives the number of blocks}
\item{param}{additional paramters specifying the distribution of the errors, the number of trials for family \code{"binomial"}}
\item{leftEnd}{a numeric vector of the same length as \code{value} containing the left end of each block}
\item{rightEnd}{a numeric vector of the same length as \code{value} containing the left end of each block}
\item{x0}{a single numeric giving the last unobserved sample point directly before sampling started, i.e. before \code{leftEnd[0]}}
\item{leftIndex}{a numeric vector of the same length as \code{value} containing the index of the sample points corresponding to the block's left end, cf. \code{\link{stepcand}}}
\item{rightIndex}{a numeric vector of the same length as \code{value} containing the index of the sample points corresponding to the block's right end, cf. \code{\link{stepcand}}}
\item{x, object}{the object}
\item{y}{a numeric vector containing the data with which to compare the fit}
\item{df}{the number of estimated parameters: by default the number of blocks for families \code{poisson} and \code{binomial}, one more (for the variance) for family \code{gauss}}
\item{nobs}{the number of observations used for estimating}
\item{\dots}{for generic methods only}
\item{i,j,drop}{see \code{"\link{[.data.frame}"}}
\item{refit}{\code{\link{logical}}; determines whether the function will be refitted after subselection, i.e. whether the selection should be interpreted as a fit with fewer jumps); in that case, for \code{family = "gaussKern"}, \code{refit} needs to be set to the original data, i.e. \code{y}}
\item{dataspace}{\code{\link{logical}} determining whether the expected value should be plotted instead of the fitted parameter value, useful e.g. for \code{family = "binomial"}, where it will plot the fitted success probability times the number of trials per observation}
}

\value{
\item{for \code{stepfit}}{an object of class \code{stepfit} which extends \code{\link{stepblock}}, additionally containing \code{\link{attr}}ibutes \code{cost}, \code{family} and \code{param}, as well as columns \code{leftIndex} and \code{rightIndex}}
\item{for \code{[.stepfit}}{an object of class \code{stepfit} which contains the selected subset}
\item{for \code{fitted.stepfit}}{a numeric vector of length \code{rightIndex[length(rightIndex)]} giving the fit at the original sample points}
\item{for \code{residuals.stepfit}}{a numeric vector of length \code{rightIndex[length(rightIndex)]} giving the residuals at the original sample points}
\item{for \code{logLik.stepfit}}{an object of class \code{\link{logLik}} giving the likelihood of the data given this fit, e.g. for use with \code{\link{AIC}} and \code{\link{stepsel}}; this will (incorrectly) treat \code{family = "gaussKern"} as if it were fitted with \code{family = "gauss"}}
\item{for \code{plot.stepfit}, \code{plot.stepfit}}{the corresponding functions for \code{\link{stepblock}} are called}
}

\seealso{\code{\link{stepblock}}, \code{\link{stepbound}}, \code{\link{steppath}}, \code{\link{stepsel}}, \link{family}, \code{"\link{[.data.frame}"}, \code{\link{fitted}}, \code{\link{residuals}}, \code{\link{logLik}}, \code{\link{AIC}}}

\examples{
# simulate 5 blocks (4 jumps) within a total of 100 data points
b <- c(sort(sample(1:99, 4)), 100)
p <- rep(runif(5), c(b[1], diff(b))) # success probabilities
# binomial observations, each with 10 trials
y <- rbinom(100, 10, p)
# find solution with 5 blocks
fit <- steppath(y, family = "binomial", param = 10)[[5]]
plot(y, ylim = c(0, 10))
lines(fit, col = "red")
# residual diagnostics for Gaussian data
yg <- rnorm(100, qnorm(p), 1)
fitg <- steppath(yg)[[5]]
plot(yg, ylim = c(0, 10))
lines(fitg, col = "red")
plot(resid(fitg, yg))
qqnorm(resid(fitg, yg))
}

\keyword{nonparametric}
