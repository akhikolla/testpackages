\name{scrHpd}
\alias{scrHpd}

\title{Calculate an SCB from a samples matrix}

\description{Calculate an SCB from a samples matrix, which minimizes
the absolute distances of the contained samples to a mode vector, at
each gridpoint. Therefore the SCB might be considered an \dQuote{HPD
SCB}.}

\usage{scrHpd(samples, mode = apply(samples, 2, median), level = 0.95)}

\arguments{
  \item{samples}{m by n matrix where m is the number of samples and n
  the number of parameters, hence each (multivariate) sample is a row in
  the matrix \code{samples}} 
  \item{mode}{mode vector of length n (defaults to the vector of medians)}
  \item{level}{credible level for the SCB (default: 0.95)}
}
\details{
  This function first computes the matrix of absolute distances of the
  samples to the mode vector. Then based on this distance matrix, a
  one-sided SCB as described in Besag et al. (1995) is computed, which
  is then mapped back to the samples.
}
\value{
  A matrix with rows \dQuote{lower} and \dQuote{upper}, with the lower
  and upper SCB bounds.
}
\references{Besag, J.; Green, P.; Higdon, D. \& Mengersen, K. (1995):
  \dQuote{Bayesian computation and stochastic systems (with
    discussion)}, \emph{Statistical Science}, 10, 3-66.}
\author{Daniel Saban\'es Bov\'e}
\seealso{\code{\link{empiricalHpd}}}
\examples{
## create some samples
time <- 1:10
nSamples <- 50
samples <- t(replicate(nSamples,
                       time * rnorm(1) + rexp(1))) +
           rnorm(length(time) * nSamples)
matplot(time, t(samples), type="l", lty=1, col=1,
        xlab="time", ylab="response")

## now test the function: 50% credible band
scb <- scrHpd(samples, level=0.5)
matlines(time, t(scb), col=2, lwd=2, lty=1)
}

\keyword{multivariate}
\keyword{htest}
