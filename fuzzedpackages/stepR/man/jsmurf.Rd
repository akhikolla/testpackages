\name{jsmurf}
\alias{jsmurf}

\title{Reconstruct filtered piecewise constant functions with noise}
\description{Reconstructs a piecewise constant function to which white noise was added and the sum filtered afterwards.

\bold{Deprecation warning:} This function is mainly used for patchlamp recordings and may be transferred to a specialised package.
}

\usage{
jsmurf(y, x = 1:length(y), x0 = 2 * x[1] - x[2], q, alpha = 0.05, r = 4e3,
  lengths = 2^(floor(log2(length(y))):floor(log2(max(length(param$kern) + 1,
    1 / param$param$cutoff)))), param, rm.out = FALSE,
    jumpint = confband, confband = FALSE)
}

\arguments{
\item{y}{a numeric vector containing the serial data}
\item{x}{a numeric vector of the same length as \code{y} containing the corresponding sample points}
\item{x0}{a single numeric giving the last unobserved sample point directly before sampling started}
\item{q}{threshold value, by default chosen automatically}
\item{alpha}{significance level; if set to a value in (0,1), \code{q} is chosen as the corresponding quantile of the asymptotic (if \code{r} is not given) null distribution (and any value specified for \code{q} is silently ignored)}
\item{r}{numer of simulations; if specified along \code{alpha}, \code{q} is chosen as the corresponding quantile of the simulated null distribution}
\item{lengths}{length of intervals considered; by default up to a sample size of 1000 all lengths, otherwise only dyadic lengths}
\item{param}{a \code{\link{dfilter}} object specifiying the filter}
\item{rm.out}{a \code{\link{logical}} specifying whether outliers should be removed prior to the analysis}
\item{jumpint}{\code{\link{logical}} (\code{FALSE} by default), indicates if confidence sets for jumps should be computed}
\item{confband}{\code{\link{logical}}, indicates if a confidence band for the piecewise-continuous function should be computed}
}

\value{
An object object of class \code{\link{stepfit}} that contains the fit; if \code{jumpint == TRUE} function \code{\link{jumpint}} allows to extract the \code{1 - alpha} confidence interval for the jumps, if \code{confband == TRUE} function \code{\link{confband}} allows to extract the \code{1 - alpha} confidence band.
}

\seealso{\code{\link{stepbound}}, \code{\link{bounds}}, \link{family}, \code{\link{MRC.asymptotic}}, \code{\link{sdrobnorm}}, \code{\link{stepfit}}}

\encoding{utf8}
\references{
Hotz, T., \enc{Schütte}{Schuette}, O., Sieling, H., Polupanow, T., Diederichsen, U., Steinem, C., and Munk, A. (2013) Idealizing ion channel recordings by a jump segmentation multiresolution filter. \emph{IEEE Transactions on NanoBioscience} \bold{12}(4), 376--386.
}

\examples{
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
sim <- contMC(0.5 * sampling, 0:1, rates, sampling=sampling, family="gaussKern",
  param = list(df=df.over, over=over, sd=0.3))
plot(sim$data, pch = ".")
lines(sim$discr, col = "red")
# fit using filter corresponding to sample rate
df <- dfilter("bessel", list(pole=4, cutoff=cutoff / sampling))
fit <- jsmurf(sim$data$y, sim$data$x, param=df, r=1e2)
lines(fit, col = "blue")
# fitted values take filter into account
lines(sim$data$x, fitted(fit), col = "green3", lty = 2)
}

\keyword{nonparametric}
