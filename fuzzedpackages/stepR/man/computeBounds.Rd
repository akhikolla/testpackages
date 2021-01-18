\name{computeBounds}
\alias{computeBounds}

\title{Computation of the bounds}
\description{Computes the multiscale contraint given by the multiscale test, (3.12) in the vignette. In more detail, returns the bounds of the interval of parameters for which the test statistic is smaller than or equal to the critical value for the corresponding length, i.e. the two solutions resulting from equating the test statistic to the critical value.\cr
If \code{q == NULL} a Monte-Carlo simulation is required for computing critical values. Since a Monte-Carlo simulation lasts potentially much longer (up to several hours or days if the number of observations is in the millions) than the main calculations, this package saves them by default in the workspace and on the file system such that a second call requiring the same Monte-Carlo simulation will be much faster. For more details, in particular to which arguments the Monte-Carlo simulations are specific, see Section \cite{Storing of Monte-Carlo simulations} below. Progress of a Monte-Carlo simulation can be reported by the argument \code{messages} and the saving can be controlled by the argument \code{option}, both can be specified in \code{\dots} and are explained in \code{\link{monteCarloSimulation}} and \code{\link{critVal}}, respectively.
}

\usage{
computeBounds(y, q = NULL, alpha = NULL, family = NULL,
              intervalSystem = NULL, lengths = NULL, ...)
}

\arguments{
\item{y}{a numeric vector containing the observations}
\item{q}{either \code{NULL}, then the vector of critical values at level \code{alpha} will be computed from a Monte-Carlo simulation, or a numeric giving the global quantile or a numeric vector giving the vector of critical values. Either \code{q} or \code{alpha} must be given. Otherwise, \code{alpha == 0.5} is chosen with a \code{\link{warning}}. This argument will be passed to \code{\link{critVal}} to obtain the needed critical values. Additional parameters for the computation of \code{q} can be specified in \code{\dots}, for more details see the documentation of \code{\link{critVal}}. Please note that by default the Monte-Carlo simulation will be saved in the workspace and on the file system, for more details see Section \cite{Storing of Monte-Carlo simulations} below}
\item{alpha}{a probability, i.e. a single numeric between 0 and 1, giving the significance level. Its choice is a trade-off between data fit and parsimony of the estimator. In other words, this argument balances the risks of missing change-points and detecting additional artefacts. For more details on this choice see (Frick et al., 2014, section 4) and (Pein et al., 2017, section 3.4). Either \code{q} or \code{alpha} must be given. Otherwise, \code{alpha == 0.5} is chosen with a \code{\link{warning}}}
\item{family}{a string specifying the assumed parametric family, for more details see \link{parametricFamily}, currently \code{"gauss"}, \code{"hsmuce"} and \code{"mDependentPS"} are supported. By default (\code{NULL}) \code{"gauss"} is assumed}
\item{intervalSystem}{a string giving the used interval system, either \code{"all"} for all intervals, \code{"dyaLen"} for all intervals of dyadic length or \code{"dyaPar"} for the dyadic partition, for more details see \link{intervalSystem}. By default (\code{NULL}) the default interval system of the specified parametric family will be used, which one this will be is described in \code{\link{parametricFamily}}}
\item{lengths}{an integer vector giving the set of lengths, i.e. only intervals of these lengths will be considered. Note that not all lengths are possible for all interval systems and for all parametric families, see \link{intervalSystem} and \link{parametricFamily}, respectively, to see which ones are allowed. By default (\code{NULL}) all lengths that are possible for the specified \code{intervalSystem} and for the specified parametric family will be used}
\item{...}{there are two groups of further arguments:
\enumerate{
\item further parameters of the parametric family. Depending on argument \code{family} some might be required, but others might be optional, please see \link{parametricFamily} for more details,
\item further parameters that will be passed to \code{\link{critVal}}. \code{\link{critVal}} will be called automatically with the number of observations \code{n = length(y)}, the arguments \code{family}, \code{intervalSystem}, \code{lengths}, \code{q} and \code{output} set. For these arguments no user interaction is required and possible, all other arguments of \code{\link{critVal}} can be passed additionally
}
}
}

\value{
A \code{\link{data.frame}} containing two integer vectors \code{li} and \code{ri} and two numeric vectors \code{lower} and \code{upper}. For each interval in the set of intervals specified by \code{intervalSystem} and \code{lengths} \code{li} and \code{ri} give the left and right index of the interval and \code{lower} and \code{upper} give the lower and upper bounds for the parameter on the given interval.
}

\section{Storing of Monte-Carlo simulations}{
If \code{q == NULL} a Monte-Carlo simulation is required for computing critical values. Since a Monte-Carlo simulation lasts potentially much longer (up to several hours or days if the number of observations is in the millions) than the main calculations, this package offers multiple possibilities for saving and loading the simulations. Progress of a simulation can be reported by the argument \code{messages} which can be specified in \code{\dots} and is explained in the documentation of \code{\link{monteCarloSimulation}}. Each Monte-Carlo simulation is specific to the number of observations, the parametric family (including certain parameters, see \link{parametricFamily}) and the interval system, and for simulations of class \code{"MCSimulationMaximum"}, additionally, to the set of lengths and the used penalty. Monte-Carlo simulations can also be performed for a (slightly) larger number of observations \eqn{n_q}{nq} given in the argument \code{nq} in \code{\dots} and explained in the documentation of \code{\link{critVal}}, which avoids extensive resimulations for only a little bit varying number of observations. Simulations can either be saved in the workspace in the variable \code{critValStepRTab} or persistently on the file system for which the package \code{\link{R.cache}} is used. Moreover, storing in and loading from variables and \link[=readRDS]{RDS} files is supported. Finally, a pre-simulated collection of simulations can be accessed by installing the package \code{stepRdata} available from \url{http://www.stochastik.math.uni-goettingen.de/stepRdata_1.0-0.tar.gz}. The simulation, saving and loading can be controlled by the argument \code{option} which can be specified in \code{\dots} and is explained in the documentation of \code{\link{critVal}}. By default simulations will be saved in the workspace and on the file system.  For more details and for how simulation can be removed see Section \cite{Simulating, saving and loading of Monte-Carlo simulations} in \code{\link{critVal}}.
}

\note{
Depending on \code{intervalSystem} and \code{lengths} the intervals might be ordered differently to allow fast computation. For most applications the order should not matter. Otherwise, the entries can be reordered with \code{\link{order}}, an example is given below.
}

\seealso{\code{\link{critVal}}, \code{\link{penalty}}, \code{\link{parametricFamily}}, \code{\link{intervalSystem}}, \code{\link{stepFit}}, \code{\link{computeStat}}, \code{\link{monteCarloSimulation}}}

\references{
Frick, K., Munk, A., Sieling, H. (2014) Multiscale change-point inference. With discussion and rejoinder by the authors. \emph{Journal of the Royal Statistical Society, Series B} \bold{76}(3), 495--580.

Pein, F., Sieling, H., Munk, A. (2017) Heterogeneous change point inference. \emph{Journal of the Royal Statistical Society, Series B}, \bold{79}(4), 1207--1227.
} 

\examples{
y <- c(rnorm(50), rnorm(50, 1))

# the multiscale contraint
bounds <- computeBounds(y, alpha = 0.5)

# the order of the bounds depends on intervalSystem and lengths
# to allow fast computation
# if a specific order is required it can be reordered by order
# b is ordered with increasing left indices and increasing right indices
b <- bounds[order(bounds$li, bounds$ri), ]
attr(b, "row.names") <- seq(along = b$li)

# higher significance level for larger detection power, but less confidence
computeBounds(y, alpha = 0.99)

# smaller significance level for stronger confidence statements, but at
# the risk of missing change-points
computeBounds(y, alpha = 0.05)
\donttest{
# different interval system, lengths, penalty and given parameter sd
computeBounds(y, alpha = 0.5, intervalSystem = "dyaLen",
              lengths = c(1L, 2L, 4L, 8L), penalty = "weights",
              weights = c(0.4, 0.3, 0.2, 0.1), sd = 0.5)
}              
# with given q
identical(computeBounds(y, q = critVal(100L, alpha = 0.5)), bounds)
identical(computeBounds(y, q = critVal(100L, alpha = 0.5, output = "value")),
          bounds)
\donttest{
# the above calls saved and (attempted to) load Monte-Carlo simulations and
# simulated them for nq = 128 observations
# in the following call no saving, no loading and simulation for n = 100
# observations is required, progress of the simulation will be reported
computeBounds(y, alpha = 0.5, messages = 1000L,
              options = list(simulation = "vector",
                             load = list(), save = list()))
                             
# with given stat to compute q
stat <- monteCarloSimulation(n = 128L)
identical(computeBounds(y, alpha = 0.5, stat = stat),
          computeBounds(y, alpha = 0.5, options = list(load = list())))
}}

\keyword{nonparametric}
