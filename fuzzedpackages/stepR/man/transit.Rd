\name{transit}
\alias{transit}

\title{TRANSIT algorithm for detecting jumps}

\description{
Reimplementation of VanDongen's algorithm for detecting jumps in ion channel recordings.

\bold{Deprecation warning:} This function is mainly used for patchlamp recordings and may be transferred to a specialised package.
}

\usage{
transit(y, x = 1:length(y), x0 = 2 * x[1] - x[2], sigma.amp = NA, sigma.slope = NA,
  amp.thresh = 3, slope.thresh = 2, rel.amp.n = 3, rel.amp.thresh = 4,
  family = c("gauss", "gaussKern"), param = NULL, refit = FALSE)
}

\arguments{
\item{y}{a numeric vector containing the serial data}
\item{sigma.amp}{amplitude (i.e. raw data within block) standard deviation; estimated using \code{\link{sdrobnorm}} if omitted}
\item{sigma.slope}{slope (i.e. central difference within block) standard deviation; estimated using \code{\link{sdrobnorm}} if omitted}
\item{amp.thresh}{amplitude threshold}
\item{slope.thresh}{slope threshold}
\item{rel.amp.n}{relative amplitude threshold will be used for blocks with no more datapoints than this}
\item{rel.amp.thresh}{relative amplitude threshold}
\item{x}{a numeric vector of the same length as \code{y} containing the corresponding sample points}
\item{x0}{a single numeric giving the last unobserved sample point directly before sampling started}
\item{family, param}{specifies distribution of data, see \link{family}}
\item{refit}{should the \code{value}s for \code{family = "gaussKern"} be obtained by fitting in the end (otherwise they are meaningless)}
}

\value{
Returns an object of class \code{\link{stepfit}} which encodes the jumps and corresponding mean values.
}

\note{
Only central, no forward differences have been used in this implementation. Moreover, the standard deviations will be estimated by \code{\link{sdrobnorm}} if omitted (respecting the filter's effect if applicable).
}

\references{
VanDongen, A. M. J. (1996) A new algorithm for idealizing single ion channel data containing multiple unknown conductance levels. \emph{Biophysical Journal} \bold{70}(3), 1303--1315.
}

\seealso{\code{\link{stepfit}}, \code{\link{sdrobnorm}}, \code{\link{jsmurf}}, \code{\link{stepbound}}, \code{\link{steppath}}}

\examples{
# estimating step-functions with Gaussian white noise added
# simulate a Gaussian hidden Markov model of length 1000 with 2 states
# with identical transition rates 0.01, and signal-to-noise ratio 2
sim <- contMC(1e3, 0:1, matrix(c(0, 0.01, 0.01, 0), 2), param=1/2)
plot(sim$data, cex = 0.1)
lines(sim$cont, col="red")
# maximum-likelihood estimation under multiresolution constraints
fit.MRC <- smuceR(sim$data$y, sim$data$x)
lines(fit.MRC, col="blue")
# choose number of jumps using BIC
path <- steppath(sim$data$y, sim$data$x, max.blocks=1e2)
fit.BIC <- path[[stepsel.BIC(path)]]
lines(fit.BIC, col="green3", lty = 2)

# estimate after filtering
# simulate filtered ion channel recording with two states
set.seed(9)
# sampling rate 10 kHz
sampling <- 1e4
# tenfold oversampling
over <- 10
# 1 kHz 4-pole Bessel-filter, adjusted for oversampling
cutoff <- 1e3
df.over <- dfilter("bessel", list(pole=4, cutoff=cutoff / sampling / over))
# two states, leaving state 1 at 10 Hz, state 2 at 20 Hz
rates <- rbind(c(0, 10), c(20, 0))
# simulate 0.5 s, level 0 corresponds to state 1, level 1 to state 2
# noise level is 0.3 after filtering
Sim <- contMC(0.5 * sampling, 0:1, rates, sampling=sampling, family="gaussKern",
  param = list(df=df.over, over=over, sd=0.3))
plot(Sim$data, pch = ".")
lines(Sim$discr, col = "red")
# fit under multiresolution constraints using filter corresponding to sample rate
df <- dfilter("bessel", list(pole=4, cutoff=cutoff / sampling))
Fit.MRC <- jsmurf(Sim$data$y, Sim$data$x, param=df, r=1e2)
lines(Fit.MRC, col = "blue")
# fit using TRANSIT
Fit.trans <- transit(Sim$data$y, Sim$data$x)
lines(Fit.trans, col = "green3", lty=2)
}

\keyword{nonparametric}
