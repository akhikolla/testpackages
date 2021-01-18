\name{simData2s}
\alias{simData2s}
\title{Simulate faecal egg count data (2-sample situation)}
\usage{
simData2s(n = 10, preMean = 500, delta = 0.1, kappa = 0.5, 
  deltaShape = NULL, phiPre = 1, phiPost = phiPre, f = 50, 
  paired = TRUE, rounding = TRUE, seed = NULL)
}
\arguments{
\item{n}{positive integer. Sample size.}

\item{preMean}{numeric. True pre-treatment epg.}

\item{delta}{numeric. Proportion of epg left after treatment, between 0 and 1. 1 - \eqn{\delta} is reduction in mean after treatment, \code{delta} = 0.1 indicates a 90\% reduction.}

\item{kappa}{numeric. Overdispersion parameter, \eqn{\kappa \to \infty}{\kappa -> \infty} corresponds to Poisson distribution.}

\item{deltaShape}{numeric. Shape parameter for the distribution of reductions. If NULL, the same reduction is applied to the latent true epg of each animal.}

\item{phiPre}{numeric. Pre-treatment prevalence (i.e. proportion of infected animals), between 0 and 1.}

\item{phiPost}{numeric. Post-treatment prevalence, between 0 and 1.}

\item{f}{integer or vector of integers. Correction factor of the egg counting technique}

\item{paired}{logical. If TRUE, paired samples are simulated. Otherwise unpaired samples are simulated.}

\item{rounding}{logical. If TRUE, the Poisson mean for the raw counts is rounded. The rounding applies since the mean epg is frequently reported as an integer value. For more information, see Details.}

\item{seed}{an integer that will be used in a call to set.seed before simulation. If NULL, a random seed is allocated.}
}
\value{
A data.frame with six columns, namely the observed epg (\code{obs}),
    actual number of eggs counted (\code{master}) and true epg in the sample (\code{true}) for both pre- and post- treatment.
}
\description{
Generates two samples of (zero-inflated) egg count data
}
\details{
In the simulation of raw (\code{master}) counts, it follows a Poisson distribution with some mean. The mean is frequently rounded down if it has a very low value and \code{rounding = TRUE}, there expects to be some bias in the mean reduction when \eqn{\mu} < 150 and \eqn{\delta} < 0.1. Set \code{rounding = FALSE} for not to have any bias.
}
\author{
Craig Wang \cr
Michaela Paul 
}
\examples{
fec <- simData2s(n = 10, preMean = 500, delta = 0.1, kappa = 0.5)

## show the bias when the true reduction should be 95\%
fec <- simData2s(n = 1e5, preMean = 150, delta = 0.05, 
          kappa = 0.5, seed = 1)
1 - mean(fec$masterPost)/mean(fec$masterPre)
## without bias
fec <- simData2s(n = 1e5, preMean = 150, delta = 0.05, 
          kappa = 0.5, seed = 1, rounding = FALSE)
1 - mean(fec$masterPost)/mean(fec$masterPre)
}
\seealso{
\code{\link{fecr_stan}} for analyzing faecal egg count data with two samples 
}
\keyword{datagen}
