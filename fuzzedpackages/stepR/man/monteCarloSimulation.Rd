\name{monteCarloSimulation}
\alias{monteCarloSimulation}

\title{Monte Carlo simulation}
\description{Performs Monte-Carlo simulations of the multiscale vector of statistics, (3.9) in the vignette, and of the penalised multiscale statistic, (3.6) in the vignette, when no signal is present, see also section 3.2.3 in the vignette.}

\usage{
monteCarloSimulation(n, r = 1e4L, family = NULL, intervalSystem = NULL,
                     lengths = NULL, penalty = NULL,
                     output = c("vector", "maximum"), seed = n,
                     rand.gen = NULL, messages = NULL, ...)
}

\arguments{
\item{n}{a positive integer giving the number of observations for which the Monte-Carlo simulation will be performed}
\item{r}{a positive integer giving the number of repititions}
\item{family}{a string specifying the assumed parametric family, for more details see \link{parametricFamily}, currently \code{"gauss"}, \code{"hsmuce"} and \code{"mDependentPS"} are supported. By default (\code{NULL}) \code{"gauss"} is assumed}
\item{intervalSystem}{a string giving the used interval system, either \code{"all"} for all intervals, \code{"dyaLen"} for all intervals of dyadic length or \code{"dyaPar"} for the dyadic partition, for more details see \link{intervalSystem}. By default (\code{NULL}) the default interval system of the specified parametric family will be used, which one this will be is described in \code{\link{parametricFamily}}}
\item{lengths}{an integer vector giving the set of lengths, i.e. only intervals of these lengths will be considered. Only required for \code{output == "maximum"}, otherwise ignored with a \code{\link{warning}}. Note that not all lengths are possible for all interval systems and for all parametric families, see \link{intervalSystem} and \link{parametricFamily}, respectively, to see which ones are allowed. By default (\code{NULL}) all lengths that are possible for the specified \code{intervalSystem} and for the specified parametric family will be used}
\item{penalty}{a string specifying how the statistics will be penalised, either \code{"sqrt"}, \code{"log"} or \code{"none"}, see \link{penalty} and section 3.2 in the vignette for more details. Only required for \code{output == "maximum"}, otherwise ignored with a \code{\link{warning}}. By default (\code{NULL}) the default penalty of the specified parametric family will be used, which one this will be is described in \code{\link{parametricFamily}}}
\item{output}{a string specifying the output, see \cite{Value}}
\item{seed}{will be passed to \code{\link{set.seed}} to set a seed, \code{\link{set.seed}} will not be called if this argument is set to \code{"no"}, i.e. a single value, interpreted as an \code{\link{integer}}, \code{NULL} or \code{"no"}}
\item{rand.gen}{by default (\code{NULL}) this argument will be replaced by the default function to generate random observations of the given \code{family}. Note that a Monte-Carlo simulation can only be saved if \code{rand.gen == NULL}. Alternatively, an own function expecting a single argument named \code{data} and returning a numeric vector of length \code{n}, this is given by \code{data$n}. Will be called with \code{rand.gen(data = data)}, with \code{data} a list containing the named entries \code{n}, the expected number of data points, and parameters of the parametric family, e.g. \code{sd} for \code{family == "gauss"} or \code{covariances} for \code{family == "mDependentPS"}}
\item{messages}{a positive integer or \code{NULL}, in each \code{messages} iteration a message will be printed in order to show the progress of the simulation, if \code{NULL} no message will be given}
\item{...}{further parameters of the parametric family. Depending on the argument \code{family} some might be required, but others might be optional, please see \link{parametricFamily} for more details}
}

\value{
If \code{output == "vector"} an object of class \code{"MCSimulationVector"}, i.e. a \eqn{d_n}{dn} times \code{r} matrix containing \code{r} independent samples of the multiscale vector of statistics, with \eqn{d_n}{dn} the number of scales, i.e. the number of possible lengths for the given interval system and given parametric family. If \code{output == "maximum"} an object of class \code{"MCSimulationMaximum"}, i.e. a vector of length \code{r} containing \code{r} independent samples of the penalised multiscale statistic. For both, additionally, the following \code{\link{attributes}} are set:
\itemize{
\item "keyList": A list specifying for which number of observations \code{n}, which parametric family with which parameters by a \acronym{SHA-1} hash, which interval system and in case of \code{"MCSimulationMaximum"}, additionally, for which lengths and which penalisation the simulation was performed.
\item "key": A key used internally for identification when the object will be saved and loaded.
\item "n": The number of observations \code{n} for which the simulation was performed.
\item "lengths": The lengths for which the simulation was performed.
\item "save": A \code{logical} which is \code{TRUE} if the object can be saved which is the case for \code{rand.gen == NULL} and \code{FALSE} otherwise. 
}
}

\seealso{\code{\link{critVal}}, \code{\link{computeStat}}, \code{\link{penalty}}, \code{\link{parametricFamily}}, \code{\link{intervalSystem}}}

\references{
Frick, K., Munk, A., Sieling, H. (2014) Multiscale change-point inference. With discussion and rejoinder by the authors. \emph{Journal of the Royal Statistical Society, Series B} \bold{76}(3), 495--580.

Pein, F., Sieling, H., Munk, A. (2017) Heterogeneous change point inference. \emph{Journal of the Royal Statistical Society, Series B}, \bold{79}(4), 1207--1227.
}

\examples{
# monteCarloSimulation will be called in critVal, can be called explicitly
# object of class MCSimulationVector
stat <- monteCarloSimulation(n = 100L)
\donttest{identical(critVal(n = 100L, alpha = 0.5, stat = stat),
          critVal(n = 100L, alpha = 0.5,
                  options = list(load = list(), simulation = "matrix")))

# object of class MCSimulationMaximum
stat <- monteCarloSimulation(n = 100L, output = "maximum")
identical(critVal(n = 100L, alpha = 0.5, stat = stat),
          critVal(n = 100L, alpha = 0.5,
                  options = list(load = list(), simulation = "vector")))

# different interval system, lengths and penalty
monteCarloSimulation(n = 100L, output = "maximum", intervalSystem = "dyaLen",
              lengths = c(1L, 2L, 4L, 8L), penalty = "log")

# with a different number of iterations, different seed,
# reported progress and user written rand.gen function
stat <- monteCarloSimulation(n = 100L, r = 1e3, seed = 1, messages = 100,
                             rand.gen = function(data) {rnorm(100)})

# the optional argument sd of parametric family "gauss" will be replaced by 1
identical(monteCarloSimulation(n = 100L, r = 1e3, sd = 5),
          monteCarloSimulation(n = 100L, r = 1e3, sd = 1))

# simulation for family "hsmuce"
monteCarloSimulation(n = 100L, family = "hsmuce")

# simulation for family "mDependentGauss"
# covariances must be given (can also be given by correlations or filter)
stat <- monteCarloSimulation(n = 100L, family = "mDependentPS",
                             covariances = c(1, 0.5, 0.3))

# variance will be standardized to 1
# output might be on some systems even identical
all.equal(monteCarloSimulation(n = 100L, family = "mDependentPS",
                               covariances = c(2, 1, 0.6)), stat)
}}

\keyword{nonparametric}
