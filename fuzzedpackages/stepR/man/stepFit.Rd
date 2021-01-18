\name{stepFit}
\alias{stepFit}

\title{Piecewise constant multiscale inference}
\description{Computes the multiscale regression estimator, see (3.1) in the vignette, and allows for confidence statements, see section 3 in the vignette. It implements the estimators \acronym{SMUCE} and \acronym{HSMUCE} as well as their confidence intervals and bands.\cr
If \code{q == NULL} a Monte-Carlo simulation is required for computing critical values. Since a Monte-Carlo simulation lasts potentially much longer (up to several hours or days if the number of observations is in the millions) than the main calculations, this package saves them by default in the workspace and on the file system such that a second call requiring the same Monte-Carlo simulation will be much faster. For more details, in particular to which arguments the Monte-Carlo simulations are specific, see Section \cite{Storing of Monte-Carlo simulations} below. Progress of a Monte-Carlo simulation can be reported by the argument \code{messages} and the saving can be controlled by the argument \code{option}, both can be specified in \code{\dots} and are explained in \code{\link{monteCarloSimulation}} and \code{\link{critVal}}, respectively.
}

\usage{
stepFit(y, q = NULL, alpha = NULL, x = 1:length(y), x0 = 2 * x[1] - x[2],
        family = NULL, intervalSystem = NULL, lengths = NULL, confband = FALSE,
        jumpint = confband, ...)
}

\arguments{
\item{y}{a numeric vector containing the observations}
\item{q}{either \code{NULL}, then the vector of critical values at level \code{alpha} will be computed from a Monte-Carlo simulation, or a numeric giving the global quantile or a numeric vector giving the vector of critical values. Either \code{q} or \code{alpha} must be given. Otherwise, \code{alpha == 0.5} is chosen with a \code{\link{warning}}. This argument will be passed to \code{\link{critVal}} to obtain the needed critical values. Additional parameters for the computation of \code{q} can be specified in \code{\dots}, for more details see the documentation of \code{\link{critVal}}. Please note that by default the Monte-Carlo simulation will be saved in the workspace and on the file system, for more details see Section \cite{Storing of Monte-Carlo simulations} below}
\item{alpha}{a probability, i.e. a single numeric between 0 and 1, giving the significance level. Its choice is a trade-off between data fit and parsimony of the estimator. In other words, this argument balances the risks of missing change-points and detecting additional artefacts. For more details on this choice see (Frick et al., 2014, section 4) and (Pein et al., 2017, section 3.4). Either \code{q} or \code{alpha} must be given. Otherwise, \code{alpha == 0.5} is chosen with a \code{\link{warning}}}
\item{x}{a numeric vector of the same length as \code{y} containing the corresponding sample points}
\item{x0}{a single numeric giving the last unobserved sample point directly before sampling started}
\item{family}{a string specifying the assumed parametric family, for more details see \link{parametricFamily}, currently \code{"gauss"}, \code{"hsmuce"} and \code{"mDependentPS"} are supported. By default (\code{NULL}) \code{"gauss"} is assumed}
\item{intervalSystem}{a string giving the used interval system, either \code{"all"} for all intervals, \code{"dyaLen"} for all intervals of dyadic length or \code{"dyaPar"} for the dyadic partition, for more details see \link{intervalSystem}. By default (\code{NULL}) the default interval system of the specified parametric family will be used, which one this will be is described in \code{\link{parametricFamily}}}
\item{lengths}{an integer vector giving the set of lengths, i.e. only intervals of these lengths will be considered. Note that not all lengths are possible for all interval systems and for all parametric families, see \link{intervalSystem} and \link{parametricFamily}, respectively, to see which ones are allowed. By default (\code{NULL}) all lengths that are possible for the specified \code{intervalSystem} and for the specified parametric family will be used}
\item{confband}{single \code{\link{logical}}, indicates if a confidence band for the piecewise-continuous function should be computed}
\item{jumpint}{single \code{\link{logical}}, indicates if confidence sets for change-points should be computed}
\item{...}{there are two groups of further arguments:
\enumerate{
\item further parameters of the parametric family. Depending on argument \code{family} some might be required, but others might be optional, please see \link{parametricFamily} for more details,
\item further parameters that will be passed to \code{\link{critVal}}. \code{\link{critVal}} will be called automatically with the number of observations \code{n = length(y)}, the arguments \code{family}, \code{intervalSystem}, \code{lengths}, \code{q} and \code{output} set. For these arguments no user interaction is required and possible, all other arguments of \code{\link{critVal}} can be passed additionally
}
}
}

\value{
An object of class \code{\link{stepfit}} that contains the fit. If \code{jumpint == TRUE} function \code{\link{jumpint}} allows to extract the \code{1 - alpha} confidence interval for the jumps. If \code{confband == TRUE} function \code{\link{confband}} allows to extract the \code{1 - alpha} confidence band.
}

\section{Storing of Monte-Carlo simulations}{
If \code{q == NULL} a Monte-Carlo simulation is required for computing critical values. Since a Monte-Carlo simulation lasts potentially much longer (up to several hours or days if the number of observations is in the millions) than the main calculations, this package offers multiple possibilities for saving and loading the simulations. Progress of a simulation can be reported by the argument \code{messages} which can be specified in \code{\dots} and is explained in the documentation of \code{\link{monteCarloSimulation}}. Each Monte-Carlo simulation is specific to the number of observations, the parametric family (including certain parameters, see \link{parametricFamily}) and the interval system, and for simulations of class \code{"MCSimulationMaximum"}, additionally, to the set of lengths and the used penalty. Monte-Carlo simulations can also be performed for a (slightly) larger number of observations \eqn{n_q}{nq} given in the argument \code{nq} in \code{\dots} and explained in the documentation of \code{\link{critVal}}, which avoids extensive resimulations for only a little bit varying number of observations. Simulations can either be saved in the workspace in the variable \code{critValStepRTab} or persistently on the file system for which the package \code{\link{R.cache}} is used. Moreover, storing in and loading from variables and \link[=readRDS]{RDS} files is supported. Finally, a pre-simulated collection of simulations can be accessed by installing the package \code{stepRdata} available from \url{http://www.stochastik.math.uni-goettingen.de/stepRdata_1.0-0.tar.gz}. The simulation, saving and loading can be controlled by the argument \code{option} which can be specified in \code{\dots} and is explained in the documentation of \code{\link{critVal}}. By default simulations will be saved in the workspace and on the file system.  For more details and for how simulation can be removed see Section \cite{Simulating, saving and loading of Monte-Carlo simulations} in \code{\link{critVal}}.
}

\seealso{\code{\link{critVal}}, \code{\link{penalty}}, \code{\link{parametricFamily}}, \code{\link{intervalSystem}}, \code{\link{monteCarloSimulation}}}

\references{
Frick, K., Munk, A., Sieling, H. (2014) Multiscale change-point inference. With discussion and rejoinder by the authors. \emph{Journal of the Royal Statistical Society, Series B} \bold{76}(3), 495--580.

Pein, F., Sieling, H., Munk, A. (2017) Heterogeneous change point inference. \emph{Journal of the Royal Statistical Society, Series B}, \bold{79}(4), 1207--1227.
}

\examples{
# generate random observations
y <- c(rnorm(50), rnorm(50, 1))
x <- seq(0.01, 1, 0.01)
plot(x, y, pch = 16, col = "grey30", ylim = c(-3, 4))

# computation of SMUCE and its confidence statements
fit <- stepFit(y, x = x, alpha = 0.5, jumpint = TRUE, confband = TRUE)
lines(fit, lwd = 3, col = "red", lty = "22")

# confidence intervals for the change-point locations
points(jumpint(fit), col = "red")
# confidence band
lines(confband(fit), lty = "22", col = "darkred", lwd = 2)

# higher significance level for larger detection power, but less confidence
stepFit(y, x = x, alpha = 0.99, jumpint = TRUE, confband = TRUE)

# smaller significance level for the small risk that the number of
# change-points is overestimated with probability not more than 5%,
# but smaller detection power
stepFit(y, x = x, alpha = 0.05, jumpint = TRUE, confband = TRUE)
\donttest{
# different interval system, lengths, penalty and given parameter sd
stepFit(y, x = x, alpha = 0.5, intervalSystem = "dyaLen",
        lengths = c(1L, 2L, 4L, 8L), penalty = "weights",
        weights = c(0.4, 0.3, 0.2, 0.1), sd = 0.5,
        jumpint = TRUE, confband = TRUE)
}        
# with given q
identical(stepFit(y, x = x, q = critVal(100L, alpha = 0.5),
                  jumpint = TRUE, confband = TRUE), fit)
identical(stepFit(y, x = x, q = critVal(100L, alpha = 0.5, output = "value"),
                  jumpint = TRUE, confband = TRUE), fit)
\donttest{
# the above calls saved and (attempted to) load Monte-Carlo simulations and
# simulated them for nq = 128 observations
# in the following call no saving, no loading and simulation for n = 100
# observations is required, progress of the simulation will be reported
stepFit(y, x = x, alpha = 0.5, jumpint = TRUE, confband = TRUE,
        messages = 1000L, options = list(simulation = "vector",
                                         load = list(), save = list()))

# with given stat to compute q
stat <- monteCarloSimulation(n = 128L)
identical(stepFit(y, x = x, alpha = 0.5, stat = stat,
                  jumpint = TRUE, confband = TRUE),
          stepFit(y, x = x, alpha = 0.5, jumpint = TRUE, confband = TRUE,
                  options = list(load = list())))
}}

\keyword{nonparametric}
