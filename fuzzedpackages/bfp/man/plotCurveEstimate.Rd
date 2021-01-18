\name{plotCurveEstimate}
\alias{plotCurveEstimate}
\alias{plotCurveEstimate.BayesMfp}
\alias{plotCurveEstimate.BmaSamples}

\title{Generic function for plotting a fractional polynomial curve estimate}
\description{
  Plot a fractional polynomial curve estimate for either a single model
  or a Bayesian model average over \code{\link{BayesMfp}} objects. Optionally,
  credible intervals and / or bands can be added to the plot.
}
\usage{
plotCurveEstimate(model, termName, plevel = 0.95, slevel = plevel,
plot = TRUE, legendPos = "topleft", rug = FALSE, partialResids=TRUE, 
hpd=TRUE,..., main = NULL) 
}

\arguments{
  \item{model}{an object of class \code{\link{BayesMfp}} or \code{\link{BmaSamples}}}
  \item{termName}{string denoting an FP term, as written by the
  \code{\link[=summary.BayesMfp]{summary}} method} 
  \item{plevel}{credible level for pointwise intervals, and \code{NULL}
  means no pointwise intervals (default: 0.95)}
  \item{slevel}{credible level for simultaneous credible band (SCB),
  \code{NULL} means no SCB (defaults to \code{plevel})}
  \item{plot}{if \code{FALSE}, only return values needed to produce the
  plot, but do not plot (default is \code{TRUE}, so a plot is made)}
  \item{legendPos}{position of coefficient estimates (for \code{\link{BayesMfp}}) or
    sample size (for \code{\link{BmaSamples}}) in the plot, \code{NULL} suppresses the
    printing (default is \dQuote{topleft})}
  \item{rug}{add a rug to the plot? (default: \code{FALSE})}
  \item{partialResids}{add partial residuals to the plot? (default:
    \code{TRUE})}
  \item{hpd}{use HPD intervals (\code{TRUE}, default) or quantile-based
  (\code{FALSE}) intervals?} 
  \item{\dots}{further arguments in case of a \code{\link{BayesMfp}}
    object (see details) and arguments for plotting with
    \code{\link{matplot}}}
  \item{main}{optional main argument for the plot}
}
\details{
  Further arguments for application on a \code{\link{BayesMfp}} object:
  \describe{
    \item{grid}{vector of unscaled abscissae, default is a length
      \code{gridSize} grid over the observed range specified by
      providing the argument \code{NULL}.}
    \item{post}{list with posterior parameters of the model, which may be
      provided manually to accelerate plotting in a loop}
    \item{gridSize}{default number of grid points used when no
      \code{grid} is supplied (default is 201)}
    \item{numSim}{number of simulations for estimation of the SCB
      (default is 500)}
  }
}
\value{
  a list of various plotting information:
  \item{original}{grid on the original covariate scale}
  \item{grid}{grid on the transformed scale}
  \item{mode}{mode curve values, only for \code{\link{BayesMfp}} object}
  \item{mean}{pointwise mean curve values, only for
    \code{\link{BmaSamples}} object}
  \item{median}{pointwise median curve values, only for
  \code{\link{BmaSamples}} object} 
  \item{plower}{lower boundaries for pointwise intervals}
  \item{pupper}{upper boundaries for pointwise intervals}
  \item{slower}{lower boundaries for SCB}
  \item{supper}{upper boundaries for SCB}
  \item{obsVals}{observed values of the covariate on the original scale}
  \item{sampleSize}{sample size underlying the curve estimate, only for
  \code{\link{BmaSamples}} object} 
  \item{partialResids}{partial residuals}
  \item{transform}{vector of shift and scale parameter}  
}

\seealso{\code{\link{BayesMfp}}, \code{\link{BmaSamples}}}
\examples{
## construct a BayesMfp object
set.seed(19)
x1 <- rnorm (n=15)
x2 <- rbinom (n=15, size=20, prob=0.5) 
x3 <- rexp (n=15)
y <- rt (n=15, df=2)

test <- BayesMfp (y ~ bfp (x2, max = 4) + uc (x1 + x3), nModels = 100,
method="exhaustive")

## plot the x2 curve estimate for the 20-th best model
p1 <- plotCurveEstimate (test[20], "x2")

## look at the returned list
str(p1)

## plot the BMA curve estimate for the same covariate
testBma <- BmaSamples (test)
p2 <- plotCurveEstimate (testBma, "x2")

## look at the returned list
str(p2)

## try the new options:
plotCurveEstimate (testBma, "x2", partialResids=FALSE, hpd=FALSE)
}

\keyword{regression}

