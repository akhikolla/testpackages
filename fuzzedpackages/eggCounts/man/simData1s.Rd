\name{simData1s}
\alias{simData1s}
\title{Simulate faecal egg count data (1-sample situation)}
\usage{
simData1s(n = 10, mean = 500, kappa = 0.5, phi = 1, 
  f = 50, rounding = TRUE, seed = NULL)
}
\arguments{
\item{n}{positive integer. Sample size.}

\item{mean}{numeric. True number of eggs per gram (epg).}

\item{kappa}{numeric. Overdispersion parameter, \eqn{\kappa \to \infty}{\kappa -> \infty} corresponds to Poisson distribution.}

\item{phi}{numeric. Prevalence, i.e. proportion of infected animals, between 0 and 1.}

\item{f}{positive integer. Correction factor of the egg counting technique, either an integer or a vector of integers with length \code{n}.}

\item{rounding}{logical. If TRUE, the Poisson mean for the raw counts is rounded. The rounding applies since the mean epg is frequently reported as an integer value. For more information, see Details.}

\item{seed}{integer. Random seed.}
}
\value{
A data.frame with three columns, namely the observed epg (\code{obs}), actual number of eggs counted (\code{master}) and
    true epg in the sample (\code{true}).
}
\description{
Simulates (zero-inflated) egg count data
}
\author{
Craig Wang \cr
Michaela Paul 
}
\details{
In the simulation of raw (\code{master}) counts, it follows a Poisson distribution with some mean. The mean is frequently rounded down if it has a very low value and \code{rounding = TRUE}, hence there expects to be a bias overall when \eqn{\mu} < 150. Set \code{rounding = FALSE} for not to have any bias in the simulated counts.
}
\examples{
fec <- simData1s(n = 10, mean = 500, 
          kappa = 0.5, phi = 0.7)
}
\seealso{
\code{\link{fec_stan}} for analyzing faecal egg count data with one sample 
}
\keyword{datagen}
