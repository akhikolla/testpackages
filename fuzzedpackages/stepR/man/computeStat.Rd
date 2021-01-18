\name{computeStat}
\alias{computeStat}

\title{Computation of the multiscale statistic}
\description{Computes the multiscale vector of penalised statistics, (3.7) in the vignette, or the penalised multiscale statistic, (3.6) in the vignette, for given signal.}

\usage{
computeStat(y, signal = 0, family = NULL, intervalSystem = NULL, lengths = NULL,
            penalty = NULL, nq = length(y),
            output = c("list", "vector", "maximum"), ...) 
}

\arguments{
\item{y}{a numeric vector containing the observations}
\item{signal}{the given signal, either a single numeric for a constant function equal to the given value or an object of class \code{\link{stepfit}}. More precisely, a \code{\link{list}} containing an integer vector \code{leftIndex}, an integer vector \code{rightIndex} and a numeric vector \code{value}, all of the same length, e.g. a \code{\link{data.frame}}, specifying a step function is enough}
\item{family}{a string specifying the assumed parametric family, for more details see \link{parametricFamily}, currently \code{"gauss"}, \code{"hsmuce"} and \code{"mDependentPS"} are supported. By default (\code{NULL}) \code{"gauss"} is assumed}
\item{intervalSystem}{a string giving the used interval system, either \code{"all"} for all intervals, \code{"dyaLen"} for all intervals of dyadic length or \code{"dyaPar"} for the dyadic partition, for more details see \link{intervalSystem}. By default (\code{NULL}) the default interval system of the specified parametric family will be used, which one this will be is described in \code{\link{parametricFamily}}}
\item{lengths}{an integer vector giving the set of lengths, i.e. only intervals of these lengths will be considered. Note that not all lengths are possible for all interval systems and for all parametric families, see \link{intervalSystem} and \link{parametricFamily}, respectively, to see which ones are allowed. By default (\code{NULL}) all lengths that are possible for the specified \code{intervalSystem} and for the specified parametric family will be used}
\item{penalty}{a string specifying how the statistics will be penalised, either \code{"sqrt"}, \code{"log"} or \code{"none"}, see \link{penalty} and section 3.2 in the vignette for more details. By default (\code{NULL}) the default penalty of the specified parametric family will be used, which one this will be is described in \code{\link{parametricFamily}}}
\item{nq}{a single integer larger than or equal to \code{length(y)} giving the number of observations used in the penalty term, see \link{penalty} for more details. The possibility to use a number larger than \code{length(y)} is given for comparisons, since a (slightly) larger number can be chosen in \code{\link{critVal}} and \code{\link{monteCarloSimulation}} to avoid extensive recomputations for (slightly) varying number of observations. For more details see also the Section \cite{Simulating, saving and loading of Monte-Carlo simulations} in \code{\link{critVal}}}
\item{output}{a string specifying the output, see \cite{Value}}
\item{...}{further parameters of the parametric family. Depending on argument \code{family} some might be required, but others might be optional, please see \link{parametricFamily} for more details}
}

\value{
If \code{output == list} a list containing in \code{maximum} the penalised multiscale statistic, i.e. the maximum over all test statistics, in \code{stat} the multiscale vector of penalised statistics, i.e. a vector of length \code{lengths} giving the maximum over all tests of that length, and in \code{lengths} the vector of lengths. If \code{output == vector} a numeric vector giving the multiscale vector of penalised statistics. If \code{output == maximum} a single numeric giving the penalised multiscale statistic. \code{-Inf} is returned for lengths for which on all intervals of that length contained in the set of intervals the \code{signal} is not constant and, hence, no test statistic can be computed. This behaves similar to \code{max(numeric(0))}.
}

\seealso{\link{parametricFamily}, \link{intervalSystem}, \link{penalty}, \code{\link{monteCarloSimulation}}, \code{\link{stepFit}}, \code{\link{computeBounds}}}

\references{
Frick, K., Munk, A., Sieling, H. (2014) Multiscale change-point inference. With discussion and rejoinder by the authors. \emph{Journal of the Royal Statistical Society, Series B} \bold{76}(3), 495--580.

Pein, F., Sieling, H., Munk, A. (2017) Heterogeneous change point inference. \emph{Journal of the Royal Statistical Society, Series B}, \bold{79}(4), 1207--1227.
}

\examples{
y <- rnorm(100)
# for the default signal = 0 a signal constant 0 is assumed
identical(computeStat(y), computeStat(y,
            signal = list(leftIndex = 1L, rightIndex = 100L, value = 0)))

# different constant value
ret <- computeStat(y, signal = 1)
# penalised multiscale statistic
identical(ret$maximum, computeStat(y, signal = 1, output = "maximum"))
# multiscale vector of penalised statistics
identical(ret$stat, computeStat(y, signal = 1, output = "vector"))

y <- c(rnorm(50), rnorm(50, 1))
# true signal
computeStat(y, signal = list(leftIndex = c(1L, 51L), rightIndex = c(50L, 100L),
                             value = c(0, 1)))

# fit satisfies the multiscale contraint, i.e.
# the penalised multiscale statistic is not larger than the used global quantile 1
computeStat(y, signal = stepFit(y, q = 1), output = "maximum") <= 1

# different interval system, lengths, penalty, given parameter sd
# and computed for an increased number of observations nq
computeStat(y, signal = list(leftIndex = c(1L, 51L), rightIndex = c(50L, 100L),
                             value = c(0, 1)), nq = 128, sd = 0.5,
            intervalSystem = "dyaLen", lengths = c(1L, 2L, 4L, 8L), penalty = "none")

# family "hsmuce"
computeStat(y, signal = mean(y), family = "hsmuce")

# family "mDependentPS"
signal <- list(leftIndex = c(1L, 13L), rightIndex = c(12L, 17L), value = c(0, -1))
y <- c(rep(0, 13), rep(-1, 4)) + 
  as.numeric(arima.sim(n = 17, list(ar = c(), ma = c(0.8, 0.5, 0.3)), sd = 1))
covariances <- as.numeric(ARMAacf(ar = c(), ma = c(0.8, 0.5, 0.3), lag.max = 3))
computeStat(y, signal = signal, family = "mDependentPS", covariances = covariances)
}

\keyword{nonparametric}
