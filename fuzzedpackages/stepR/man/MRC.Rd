\name{MRC}
\alias{MRC}
\alias{MRCoeff}
\alias{MRC.pvalue}
\alias{MRC.quant}
\alias{MRC.simul}
\alias{chi}
\alias{chi.FFT}
\alias{MRC.FFT}
\alias{MRCoeff.FFT}
\alias{kMRC.pvalue}
\alias{kMRC.quant}
\alias{kMRC.simul}

\title{Compute Multiresolution Criterion}

\description{
Computes multiresolution coefficients, the corresponding criterion, simulates these for Gaussian white or coloured noise, based on which p-values and quantiles are obtained.

\bold{Deprecation warning:} The function \code{MRC.simul} is deprecated, but still working, however, may be defunct in a future version. Please use instead the function \code{\link{monteCarloSimulation}}. An example how to reproduce results is given below. Some other functions are help function and might be removed, too.
}

\usage{
MRC(x, lengths = 2^(floor(log2(length(x))):0), norm = sqrt(lengths),
  penalty = c("none", "log", "sqrt"))
MRCoeff(x, lengths = 2^(floor(log2(length(x))):0), norm = sqrt(lengths), signed = FALSE)
MRC.simul(n, r, lengths = 2^(floor(log2(n)):0), penalty = c("none", "log", "sqrt"))
MRC.pvalue(q, n, r, lengths = 2^(floor(log2(n)):0), penalty = c("none", "log", "sqrt"),
  name = ".MRC.table", pos = .GlobalEnv, inherits = TRUE)
MRC.FFT(epsFFT, testFFT, K = matrix(TRUE, nrow(testFFT), ncol(testFFT)), lengths,
  penalty = c("none", "log", "sqrt"))
MRC.quant(p, n, r, lengths = 2^(floor(log2(n)):0), penalty = c("none", "log", "sqrt"),
  name = ".MRC.table", pos = .GlobalEnv, inherits = TRUE, \dots)
kMRC.simul(n, r, kern, lengths = 2^(floor(log2(n)):ceiling(log2(length(kern)))))
kMRC.pvalue(q, n, r, kern, lengths = 2^(floor(log2(n)):ceiling(log2(length(kern)))),
  name = ".MRC.ktable", pos = .GlobalEnv, inherits = TRUE)
kMRC.quant(p, n, r, kern, lengths = 2^(floor(log2(n)):ceiling(log2(length(kern)))),
  name = ".MRC.ktable", pos = .GlobalEnv, inherits = TRUE, \dots)
}

\arguments{
\item{x}{a vector of numerical observations}
\item{lengths}{vector of interval lengths to use, dyadic intervals by default}
\item{signed}{whether signed coefficients should be returned}
\item{q}{quantile}
\item{n}{length of data set}
\item{r}{number of simulations to use}
\item{name,pos,inherits}{under which name and where precomputed results are stored, or retrieved, see \code{\link{assign}}}
\item{K}{a \code{\link{logical}} matrix indicating the set of valid intervals}
\item{epsFFT}{a vector containg the FFT of the data set}
\item{testFFT}{a matrix containing the FFTs of the intervals}
\item{kern}{a filter kernel}
\item{penalty}{penalty term in the multiresolution statistic: \code{"none"} for no penalty, \code{"log"} for penalizing the log-length of an interval, and \code{"sqrt"} for penalizing the square root of the MRC; or a \code{\link{function}} taking two arguments, the first being the multiresolution coefficients, the second the interval lenghts}
\item{norm}{how the partial sums should be normalised, by default \code{sqrt(lengths)}, so they are normalised to equal variance across all interval lengths}
\item{p}{p-value}
\item{\dots}{further arguments passed to function \code{\link{quantile}}}
}

\value{
\item{for \code{MRC}}{a vector giving the maximum as well as the indices of the corresponding interval's start and length}
\item{for \code{MRCoeff}}{a matrix giving the multiresolution coefficients for all test intervals}
\item{for \code{MRC.pvalue} / \code{MRC.quant} / \code{MRC.simul}}{the corresponding p-value / quantile / vector of simulated values under the assumption of standard Gaussian white noise}
\item{for \code{kMRC.pvalue} / \code{kMRC.simul} / \code{kMRC.simul}}{the corresponding p-value / quantile / vector of simulated values under the assumption of filtered Gaussian white noise}
}

\seealso{\code{\link{monteCarloSimulation}}, \code{\link{smuceR}}, \code{\link{jsmurf}}, \code{\link{stepbound}}, \code{\link{stepsel}}, \code{\link{quantile}}}

\encoding{utf8}
\references{
Davies, P. L., Kovac, A. (2001) Local extremes, runs, strings and multiresolution. \emph{The Annals of Statistics} \bold{29}, 1--65.

\enc{Dümbgen}{Duembgen}, L., Spokoiny, V. (2001) Multiscale testing of qualitative hypotheses. \emph{The Annals of Statistics} \bold{29}, 124--152.

Siegmund, D. O., Venkatraman, E. S. (1995) Using the generalized likelihood ratio statistic for sequential detection of a change-point. \emph{The Annals of Statistics} \bold{23}, 255--271.

Siegmund, D. O., Yakir, B. (2000) Tail probabilities for the null distribution of scanning statistics. \emph{Bernoulli} \bold{6}, 191--213.
}

\examples{
set.seed(100)
all.equal(MRC.simul(100, r = 100),
          sort(monteCarloSimulation(n = 100, r = 100, output = "maximum",
                                    penalty = "none", intervalSystem = "dyaLen")),
          check.attributes = FALSE)

# simulate signal of 100 data points
set.seed(100)
f <- rep(c(0, 2, 0), c(60, 10, 30))
# add gaussian noise
x <- f + rnorm(100)
# compute multiresolution criterion
m <- MRC(x)
# compute Monte-Carlo p-value based on 100 simulations
MRC.pvalue(m["max"], length(x), 100)
# compute multiresolution coefficients
M <- MRCoeff(x)
# plot multiresolution coefficients, colours show p-values below 5% in 1% steps
op <- par(mar = c(5, 4, 2, 4) + 0.1)
image(1:length(x), seq(min(x), max(x), length = ncol(M)), apply(M[,ncol(M):1], 1:2,
  MRC.pvalue, n = length(x), r = 100), breaks = (0:5) / 100,
  col = rgb(1, seq(0, 1, length = 5), 0, 0.75),
  xlab = "location / left end of interval", ylab ="measurement",
  main = "Multiresolution Coefficients",
  sub = paste("MRC p-value =", signif(MRC.pvalue(m["max"], length(x), 100), 3)))
axis(4, min(x) + diff(range(x)) * ( pretty(1:ncol(M) - 1) ) / dim(M)[2],
  2^pretty(1:ncol(M) - 1))
mtext("interval lengths", 4, 3)
# plot signal and its mean
points(x)
lines(f, lty = 2)
abline(h = mean(x))
par(op)
}

\keyword{nonparametric}
