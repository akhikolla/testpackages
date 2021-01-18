\name{jumpint}
\alias{jumpint}
\alias{jumpint.stepfit}
\alias{points.jumpint}
\alias{confband}
\alias{confband.stepfit}
\alias{lines.confband}

\title{Confidence intervals for jumps and confidence bands for step functions}

\description{
Extract and plot confidence intervals and bands from fits given by a \code{\link{stepfit}} object.
}

\usage{
jumpint(sb, \dots)
\method{jumpint}{stepfit}(sb, \dots)
\method{points}{jumpint}(x, pch.left = NA, pch.right = NA, y.left = NA, y.right = NA, xpd = NA, \dots)
confband(sb, \dots)
\method{confband}{stepfit}(sb, \dots)
\method{lines}{confband}(x, dataspace = TRUE, \dots)
}

\arguments{
\item{sb}{the result of a fit by \code{\link{stepbound}}}
\item{x}{the object}
\item{pch.left, pch.right}{the plotting character to use for the left/right end of the interval with defaults \code{"("} and \code{"]"} (see parameter \code{pch} of \code{\link{par}})}
\item{y.left, y.right}{at which height to plot the interval boundaries with default \code{\link{par}()$usr[3]}}
\item{xpd}{see \code{\link{par}}}
\item{dataspace}{\code{\link{logical}} determining whether the expected value should be plotted instead of the fitted parameter value, useful e.g. for \code{family = "binomial"}, where it will plot the fitted success probability times the number of trials per observation}
\item{\dots}{arguments to be passed to generic methods}
}

\value{
For \code{jumpint} an object of class \code{jumpint}, i.e. a \code{\link{data.frame}} whose columns \code{rightEndLeftBound} and \code{rightEndRightBound} specify the left and right end of the confidence interval for the block's right end, resp., given the number of blocks was estimated correctly, and similarly columns \code{rightIndexLeftBound} and \code{rightIndexRightBound} specify the left and right indices of the confidence interval, resp. Function \code{\link{points}} plots these intervals on the lower horizontal axis (by default).

For \code{confband} an object of class \code{confband}, i.e. a \code{\link{data.frame}} with columns \code{lower} and \code{upper} specifying a confidence band computed at every point \code{x}; this is a simultaneous confidence band assuming the true number of jumps has been determined. Function \code{\link{lines}} plots the confidence band.
}

\note{
Observe that jumps may occur immediately before or after an observed \code{x}; this lack of knowledge is reflected in the visual impressions by the lower and upper envelopes jumping vertically early, so that possible jumps between \code{x}s remain within the band, and by the confidence intervals starting immediately after the last \code{x} for which there cannot be a jump, cf. the note in the help for \code{\link{stepblock}}.
}

\seealso{\code{\link{stepbound}}, \code{\link{points}}, \code{\link{lines}}}

\examples{
# simulate Bernoulli data with four blocks
y <- rbinom(200, 1, rep(c(0.1, 0.7, 0.3, 0.9), each=50))
# fit step function
sb <- stepbound(y, family="binomial", param=1, confband=TRUE)
plot(y, pch="|")
lines(sb)
# confidence intervals for jumps
jumpint(sb)
points(jumpint(sb), col="blue")
# confidence band
confband(sb)
lines(confband(sb), lty=2, col="blue")
}

\keyword{nonparametric}
