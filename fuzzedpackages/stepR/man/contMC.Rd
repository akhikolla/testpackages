\name{contMC}
\alias{contMC}

\title{Continuous time Markov chain}

\description{
Simulate a continuous time Markov chain.

\bold{Deprecation warning:} This function is mainly used for patchlamp recordings and may be transferred to a specialised package.
}

\usage{
contMC(n, values, rates, start = 1, sampling = 1, family = c("gauss", "gaussKern"),
  param = NULL)
}

\arguments{
\item{n}{number of data points to simulate}
\item{values}{a \code{\link{numeric}} vector specifying signal amplitudes for different states}
\item{rates}{a square \code{\link{matrix}} matching the dimension of \code{values} each with \code{rates[i,j]} specifying the transition rate from state \code{i} to state \code{j}; the diagonal entries are ignored}
\item{start}{the state in which the Markov chain is started}
\item{sampling}{the sampling rate}
\item{family}{whether Gaussian white (\code{"gauss"}) or coloured (\code{"gaussKern"}), i.e. filtered, noise should be added; cf. \link{family}}
\item{param}{for \code{family="gauss"}, a single non-negative \code{\link{numeric}} specifying the standard deviation of the noise; for \code{family="gaussKern"}, \code{param} must be a list with entry \code{df} giving the \code{\link{dfilter}} object used for filtering, an \code{\link{integer}} entry \code{over} which specifies the oversampling factor of the filter, i.e. \code{param$df} has to be created for a sampling rate of \code{sampling} times \code{over}, and an additional non-negative \code{\link{numeric}} entry \code{sd} specifying the noise's standard deviation \emph{after} filtering; cf. \link{family}}
}

\value{
A \code{\link{list}} with components
\item{\code{cont}}{an object of class \code{\link{stepblock}} containing the simulated true values in continuous time, with an additional column \code{state} specifying the corresponding state}
\item{\code{discr}}{an object of class \code{\link{stepblock}} containing the simulated true values reduced to discrete time, i.e. containing only the observable blocks}
\item{\code{data}}{a \code{\link{data.frame}} with columns \code{x} and \code{y} containing the times and values of the simulated observations, respectively}
}

\note{
This follows the description for simulating ion channels given by VanDongen (1996).
}

\references{
VanDongen, A. M. J. (1996) A new algorithm for idealizing single ion channel data containing multiple unknown conductance levels. \emph{Biophysical Journal} \bold{70}(3), 1303--1315.
}

\seealso{\code{\link{stepblock}}, \code{\link{jsmurf}}, \code{\link{stepbound}}, \code{\link{steppath}}, \code{\link{family}}, \code{\link{dfilter}}}

\examples{
# Simulate filtered ion channel recording with two states
set.seed(9)
# sampling rate 10 kHz
sampling <- 1e4
# tenfold oversampling
over <- 10
# 1 kHz 4-pole Bessel-filter, adjusted for oversampling
cutoff <- 1e3
df <- dfilter("bessel", list(pole=4, cutoff=cutoff / sampling / over))
# two states, leaving state 1 at 1 Hz, state 2 at 10 Hz
rates <- rbind(c(0, 1e0), c(1e1, 0))
# simulate 5 s, level 0 corresponds to state 1, level 1 to state 2
# noise level is 0.1 after filtering
sim <- contMC(5 * sampling, 0:1, rates, sampling=sampling, family="gaussKern",
  param = list(df=df, over=over, sd=0.1))
sim$cont
plot(sim$data, pch = ".")
lines(sim$discr, col = "red")
# noise level after filtering, estimated from first block
sd(sim$data$y[1:sim$discr$rightIndex[1]])
# show autocovariance in first block
acf(ts(sim$data$y[1:sim$discr$rightIndex[1]], freq=sampling), type = "cov")
# power spectrum in first block
s <- spec.pgram(ts(sim$data$y[1:sim$discr$rightIndex[1]], freq=sampling), spans=c(200,90))
# cutoff frequency is where power spectrum is halved
abline(v=cutoff, h=s$spec[1] / 2, lty = 2)
}

\keyword{nonparametric}
