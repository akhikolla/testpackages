\name{critVal}
\alias{critVal}

\title{Critical values}
\description{Computes the vector of critical values or the global quantile. This function offers two ways of computation, either at significance level \code{alpha} from a Monte-Carlo simulation, see also section 3.2 in the vignette for more details, or from the global quantile / critical values given in the argument \code{q}. For more details on these two options see Section \cite{Computation of critical values / global quantile}.\cr
Since a Monte-Carlo simulation lasts potentially much longer (up to several hours or days if the number of observations is in the millions) than the main calculations, this package saves them by default in the workspace and on the file system such that a second call requiring the same Monte-Carlo simulation will be much faster. For more details, in particular to which arguments the Monte-Carlo simulations are specific, see Section \cite{Storing of Monte-Carlo simulations} below. Progress of a Monte-Carlo simulation can be reported by the argument \code{messages} in \code{\dots}, explained in \code{\link{monteCarloSimulation}}, and the saving can be controlled by the argument \code{option}.
}

\usage{
critVal(n, q = NULL, alpha = NULL, nq = 2L^(as.integer(log2(n) + 1e-12) + 1L) - 1L,
        family = NULL, intervalSystem = NULL, lengths = NULL, penalty = NULL,
        weights = NULL, stat = NULL, r = 1e4, output = c("vector", "value"),
        options = NULL, ...)
}

\arguments{
\item{n}{a positive integer giving the number of observations}
\item{q}{either \code{NULL}, then the vector of critical values at level \code{alpha} will be computed from a Monte-Carlo simulation, or a numeric giving the global quantile or a numeric vector giving the vector of critical values. For more detailed information, in particular of which length the numeric vector should be, see Section \cite{Computation of critical values / global quantile}. Either \code{q} or \code{alpha} must be given. Otherwise, \code{alpha == 0.5} is chosen with a \code{\link{warning}}. Please note that by default the Monte-Carlo simulation will be saved in the workspace and on the file system, for more details see Section \cite{Simulating, saving and loading of Monte-Carlo simulations} below}
\item{alpha}{a probability, i.e. a single numeric between 0 and 1, giving the significance level. Its choice is a trade-off between data fit and parsimony of the estimator. In other words, this argument balances the risks of missing change-points and detecting additional artefacts. For more details on this choice see (Frick et al., 2014, section 4) and (Pein et al., 2017, section 3.4). Either \code{q} or \code{alpha} must be given. Otherwise, \code{alpha == 0.5} is chosen with a \code{\link{warning}}}
\item{nq}{a positive integer larger than or equal to \code{n} giving the (increased) number of observations for the Monte-Carlo simulation. See Section \cite{Simulating, saving and loading of Monte-Carlo simulations} for more details}
\item{family}{a string specifying the assumed parametric family, for more details see \link{parametricFamily}, currently \code{"gauss"}, \code{"hsmuce"} and \code{"mDependentPS"} are supported. By default (\code{NULL}) \code{"gauss"} is assumed}
\item{intervalSystem}{a string giving the used interval system, either \code{"all"} for all intervals, \code{"dyaLen"} for all intervals of dyadic length or \code{"dyaPar"} for the dyadic partition, for more details see \link{intervalSystem}. By default (\code{NULL}) the default interval system of the specified parametric family will be used, which one this will be is described in \code{\link{parametricFamily}}}
\item{lengths}{an integer vector giving the set of lengths, i.e. only intervals of these lengths will be considered. Note that not all lengths are possible for all interval systems and for all parametric families, see \link{intervalSystem} and \link{parametricFamily}, respectively, to see which ones are allowed. By default (\code{NULL}) all lengths that are possible for the specified \code{intervalSystem} and for the specified parametric family will be used}
\item{penalty}{a string specifying how different scales will be balanced, either \code{"sqrt"}, \code{"weights"}, \code{"log"} or \code{"none"}, see \link{penalty} and section 3.2 in the vignette for more details. By default (\code{NULL}) the default penalty of the specified parametric family will be used, which one this will be is described in \code{\link{parametricFamily}}}
\item{weights}{a numeric vector of length \code{length(lengths)} with only positive entries giving the weights that will be used for penalty \code{"weights"}, see \link{penalty} and section 3.2.2 in the vignette for more details. By default (\code{NULL}) equal weights will be used, i.e. \preformatted{weights == rep(1 / length(lengths), length(lengths))}}
\item{stat}{an object of class \code{"MCSimulationVector"} or \code{"MCSimulationMaximum"} giving a Monte-Carlo simulations, usually computed by \code{\link{monteCarloSimulation}}. If \code{penalty == "weights"} only \code{"MCSimulationVector"} is allowed. Has to be simulated for at least the given number of observations \code{n} and for the given \code{family}, \code{intervalSystem} and if \code{"MCSimulationMaximum"} for the given \code{lengths} and \code{penalty}. By default (\code{NULL}) the required simulation will be made available automatically accordingly to the given \code{options}. For more details see Section \cite{Simulating, saving and loading of Monte-Carlo simulations} and section 3.4 in the vignette}
\item{r}{a positive integer giving the required number of Monte-Carlo simulations if they will be simulated or loaded from the workspace or the file system}
\item{output}{a string specifying the return value, if \code{output == "vector"} the vector of critical values will be computed and if \code{output == "value"} the global quantile will be computed. For \code{penalty == "weights"} the output must be \code{"vector"}, since no global quantile can be determined for this penalty}
\item{options}{a \code{\link{list}} specifying how Monte-Carlo simulations will be simulated, saved and loaded. For more details see Section \cite{Simulating, saving and loading of Monte-Carlo simulations} and section 3.4 in the vignette}
\item{...}{there are two groups of further arguments:
\itemize{
\item further parameters of the parametric family. Depending on the argument \code{family} some might be required, but others might be optional, please see \link{parametricFamily} for more details
\item further arguments (\code{seed}, \code{rand.gen} and \code{messages}) that will be passed to \code{\link{monteCarloSimulation}}. \code{\link{monteCarloSimulation}} will be called automatically and most of the arguments will be set accordingly to the arguments of \code{critVal}, no user interaction is required and possible for these parameters. In addition, \code{seed}, \code{rand.gen} and \code{messages} can be passed by the user
}}
}

\value{
If \code{output == "vector"} a numeric vector giving the vector of critical values, i.e. a vector of length \code{length(lengths)}, giving for each length the corresponding critical value. If \code{output == "value"} a single numeric giving the global quantile. In both cases, additionally, an \code{\link[=attributes]{attribute}} \code{"n"} gives the number of observations for which the Monte-Carlo simulation was performed.
}

\section{Computation of critical values / global quantile}{
This function offers two ways to compute the resulting value:
\itemize{
\item If \code{q == NULL} it will be computed at significance level \code{alpha} from a Monte-Carlo simulation. For penalties \code{"sqrt"}, \code{"log"} and \code{"none"} the global quantile will be the (\code{1-alpha})-quantile of the penalised multiscale statistic, see section 3.2.1 in the vignette. And if required the vector of critical values will be derived from it. For penalty \code{"weights"} the vector of critical values will be calculated accordingly to the given \code{weights}. The Monte-Carlo simulation can either be given in \code{stat} or will be attempted to load or will be simulated. How Monte-Carlo simulations are simulated, saved and loaded can be controlled by the argument \code{option}, for more details see the Section \cite{Simulating, saving and loading of Monte-Carlo simulations}.
\item If \code{q} is given it will be derived from it. For the argument \code{q} either a single finite numeric giving the global quantile or a vector of finite numerics giving the vector of critical values (not allowed for \code{output == "value"}) is possible:
\itemize{
\item A single numeric giving the global quantile. If \code{output == "vector"} the vector of critical values will be computed from it for the given \code{lengths} and \code{penalty} (penalty \code{"weights"} is not allowed). Note that the global quantile is specific to the arguments \code{family}, \code{intervalSystem}, \code{lengths} and \code{penalty}.
\item A vector of length \code{length(lengths)}, giving for each length the corresponding critical value. This vector is identical to the vector of critical values.
\item A vector of length \code{n} giving for each length \code{1:n} the corresponding critical value.
\item A vector of length equal to the number of all possible lengths for the given interval system and the given parametric family giving for each possible length the corresponding critical value.
}
Additionally, an \code{\link[=attributes]{attribute}} \code{"n"} giving the number of observations for which \code{q} was computed is allowed. This \code{\link[=attributes]{attribute}} must be a single integer and equal to or larger than the argument \code{n} which means that \code{q} must have been computed for at least \code{n} observations. This allows additionally:
\itemize{
\item A vector of length \code{attr(q, "n")} giving for each length \code{1:attr(q, "n")} the corresponding critical value.
\item A vector of length equal to the number of all possible lengths for the given interval system and the given parametric family if the number of observations is \code{attr(q, "n")} giving for each possible length the corresponding critical value.
}
The \code{\link[=attributes]{attribute}} \code{"n"} will be kept or set to \code{n} if missing.
}
}

\section{Simulating, saving and loading of Monte-Carlo simulations}{
Since a Monte-Carlo simulation lasts potentially much longer (up to several hours or days if the number of observations is in the millions) than the main calculations, this function offers multiple possibilities for saving and loading the simulations. The simulation, saving and loading can be controlled by the argument \code{option}. This argument has to be a \code{\link{list}} or \code{NULL} and the following named entries are allowed: \code{"simulation"}, \code{"save"}, \code{"load"}, \code{"envir"} and \code{"dirs"}. All missing entries will be set to their default option.\cr
Objects of class \code{"MCSimulationVector"}, containing simulations of the multiscale vector of statistics, and objects of class \code{"MCSimulationMaximum"}, containing simulations of the penalised multiscale statistic (for penalties \code{"sqrt"}, \code{"log"} and \code{"none"}), can be simulated, saved and loaded. Each Monte-Carlo simulation is specific to the number of observations, the parametric family and the interval system, for \code{"MCSimulationMaximum"} additionally to the set of lengths and the used penalty. Both types will lead to the same result, however, an object of class \code{"MCSimulationVector"} is more flexible, since critical values for all penalties and all set of lengths can be derived from it, but requires much more storage space and has slightly larger saving and loading times. Note that Monte-Carlo simulations can only be saved and loaded if they are generated with the default function for generating random observations, i.e. when \code{rand.gen} (in \code{\dots}) is \code{NULL}. For a given simulation this is signalled by the \code{\link[=attributes]{attribute}} \code{"save"} which is \code{TRUE} if a simulation can be saved and \code{FALSE} otherwise.\cr
Monte-Carlo simulations can also be performed for a (slightly) larger number of observations \eqn{n_q}{nq} given in the argument \code{nq}, which avoids extensive resimulations for only a little bit varying number of observations. The overestimation control is still satisfied but the detection power is (slightly) smaller. But note that the default \code{lengths} might change when the number of observations is increased and, hence, for type \code{"vectorIncreased"} still a different simulation might be required.\cr
We refer to the different types as follow:
\itemize{
\item \code{"vector"}: an object of class \code{"MCSimulationMaximum"}, i.e. simulations of the penalized multiscale statistic, for \code{n} observations
\item \code{"vectorIncreased"}: an object of class \code{"MCSimulationMaximum"}, i.e. simulations of the penalized multiscale statistic, for \code{nq} observations
\item \code{"matrix"}: an object of class \code{"MCSimulationVector"}, i.e. simulations of the multiscale vector of statistics, for \code{n} observations
\item \code{"matrixIncreased"}: an object of class \code{"MCSimulationVector"}, i.e. simulations of the multiscale vector of statistics, for \code{nq} observations
}
The simulations can either be saved in the workspace in the variable \code{critValStepRTab} or persistently on the file system for which the package \code{\link{R.cache}} is used. Loading from the workspace is faster, but either the user has to store the workspace manually or in a new session simulations have to be performed again. Moreover, storing in and loading from variables and \link[=readRDS]{RDS} files is supported. Finally, a pre-computed collection of simulations of type \code{"matrixIncreased"} for parametric families \code{"gauss"} and \code{"hsmuce"} can be accessed by installing the package \code{stepRdata} available from \url{http://www.stochastik.math.uni-goettingen.de/stepRdata_1.0-0.tar.gz}.\cr

\subsection{options$envir and options$dirs}{
For loading from / saving in the workspace the variable \code{critValStepRTab} in the \code{\link{environment}} \code{options$envir} will be looked for and if missing in case of saving also created there. Moreover, the variable(s) specified in \code{options$save$variable} (explained in the Subsection \cite{Saving: options$save}) will be assigned to this \code{\link{environment}}. By default, the \link[=.GlobalEnv]{global environment} \code{.GlobalEnv} is used, i.e. \code{options$envir == .GlobalEnv}.\cr
For loading from / saving on the file system \code{\link{loadCache}(key = keyList, dirs = options$dirs)} and \code{\link{saveCache}(stat, key = attr(stat, "keyList"), dirs = options$dirs)} are called, respectively. In other words, \code{options$dirs} has to be a \code{\link{character}} \code{\link{vector}} constituting the path to the cache subdirectory relative to the cache root directory as returned by \code{\link{getCacheRootPath}}(). If \code{options$dirs == ""} the path will be the cache root path. By default the subdirectory \code{"stepR"} is used, i.e. \code{options$dirs == "stepR"}. Missing directories will be created.
}

\subsection{Simulation: options$simulation}{
Whenever Monte-Carlo simulations have to be performed, i.e. when \code{stat == NULL} and the required Monte-Carlo simulation could not be loaded, the type specified in \code{options$simulation} will be simulated by \code{\link{monteCarloSimulation}}. In other words, \code{options$simulation} must be a single string of the following: \code{"vector"}, \code{"vectorIncreased"}, \code{"matrix"} or \code{"matrixIncreased"}. By default (\code{options$simulation == NULL}), an object of class \code{"MCSimulationVector"} for \code{nq} observations will be simulated, i.e. \code{options$simulation} \code{== "matrixIncreased"}. For this choice please recall the explanations regarding computation time and flexibility at the beginning of this section.
}

\subsection{Loading: options$load}{
Loading of the simulations can be controlled by the entry \code{options$load} which itself has to be a \code{\link{list}} with possible entries: \code{"RDSfile"}, \code{"workspace"}, \code{"package"} and \code{"fileSystem"}. Missing entries disable the loading from this option.
Whenever a Monte-Carlo simulation is required, i.e. when the variable \code{q} is not given, it will be searched for at the following places in the given order until found:
\enumerate{
\item in the variable \code{stat},
\item in \code{options$load$RDSfile} as an \link[=readRDS]{RDS} file, i.e. the simulation will be loaded by \preformatted{readRDS(options$load$RDSfile).} In other words, \code{options$load$RDSfile} has to be a \code{\link{connection}} or the name of the file where the \R object is read from,
\item in the workspace or on the file system in the following order: \code{"vector"}, \code{"matrix"}, \code{"vectorIncreased"} and finally of \code{"matrixIncreased"}. For \code{penalty == "weights"} it will only be looked for \code{"matrix"} and \code{"matrixIncreased"}. For each options it will first be looked in the workspace and then on the file system. All searches can be disabled by not specifying the corresponding string in \code{options$load$workspace} and \code{options$load$fileSystem}. In other words, \code{options$load$workspace} and \code{options$load$fileSystem} have to be vectors of strings containing none, some or all of \code{"vector"}, \code{"matrix"}, \code{"vectorIncreased"} and \code{"matrixIncreased"},
\item in the package \code{stepRdata} (if installed) and if \code{options$load$package == TRUE}. In other words, \code{options$load$package} must be a single logical or \code{NULL},
\item if all other options fail a Monte-Carlo simulation will be performed.
}
By default (if \code{options$load} is missing / \code{NULL}) no \link[=readRDS]{RDS} file is specified and all other options are enabled, i.e. \preformatted{
options$load <- list(workspace = c("vector", "vectorIncreased",
                                   "matrix", "matrixIncreased"),
                     fileSystem = c("vector", "vectorIncreased",
                                    "matrix", "matrixIncreased"),
                     package = TRUE, RDSfile = NULL).}
}

\subsection{Saving: options$save}{
Saving of the simulations can be controlled by the entry \code{options$save} which itself has to be a \code{\link{list}} with possible entries: \code{"workspace"}, \code{"fileSystem"}, \code{"RDSfile"} and \code{"variable"}. Missing entries disable the saving in this option.\cr
All available simulations, no matter whether they are given by \code{stat}, loaded, simulated or in case of \code{"vector"} and \code{"vectorIncreased"} computed from \code{"matrix"} and \code{"matrixIncreased"}, respectively, will be saved in all options for which the corresponding type is specified. Here we say a simulation is of type \code{"vectorIncreased"} or \code{"matrixIncreased"} if the simulation is not performed for \code{n} observations. More specifically, a simulation will be saved:
\enumerate{
\item in the workspace or on the file system if the corresponding string is contained in \code{options$save$workspace} and \code{options$save$fileSystem}, respectively. In other words, \code{options$save$workspace} and \code{options$save$fileSystem} have to be vectors of strings containing none, some or all of \code{"vector"}, \code{"matrix"}, \code{"vectorIncreased"} and \code{"matrixIncreased"},
\item in an \link[=readRDS]{RDS} file specified by \code{options$save$RDSfile} which has to be a vector of one or two \code{\link[=connection]{connections}} or names of files where the \R object is saved to. If \code{options$save$RDSfile} is of length two a simulation of type \code{"vector"} or \code{"vectorIncreased"} (only one can occur at one function call) will be saved in \code{options$save$RDSfile[1]} by \preformatted{saveRDS(stat, file = options$save$RDSfile[1])} and \code{"matrix"} or \code{"matrixIncreased"} (only one can occur at one function call) will be saved in \code{options$save$RDSfile[2]}. If \code{options$save$RDSfile} is of length one both will be saved in \code{options$save$RDSfile} which means if both occur at the same call only \code{"vector"} or \code{"vectorIncreased"} will be saved. Each saving can be disabled by not specifying \code{options$save$RDSfile} or by passing an empty string to the corresponding entry of \code{options$save$RDSfile}.
\item in a variable named by \code{options$save$variable} in the \code{\link{environment}} \code{options$envir}. Hence, \code{options$save$variable} has to be a vector of one or two containing variable names (character vectors). If \code{options$save$variable} is of length two a simulation of type \code{"vector"} or \code{"vectorIncreased"} (only one can occur at one function call) will be saved in \code{options$save$variable[1]} and \code{"matrix"} or \code{"matrixIncreased"} (only one can occur at one function call) will be saved in \code{options$save$variable[2]}. If \code{options$save$variable} is of length one both will be saved in \code{options$save$variable} which means if both occur at the same call only \code{"vector"} or \code{"vectorIncreased"} will be saved. Each saving can be disabled by not specifying \code{options$save$variable} or by passing \code{""} to the corresponding entry of \code{options$save$variable}.
}
By default (if \code{options$save} is missing) \code{"vector"} and \code{"vectorIncreased"} will be saved in the workspace and \code{"matrix"} and \code{"matrixIncreased"} on the file system, i.e. \preformatted{
options$save <- list(workspace = c("vector", "vectorIncreased"),
                     fileSystem = c("matrix", "matrixIncreased"),
                     RDSfile = NULL, variable = NULL).}
Simulations can be removed from the workspace by removing the variable \code{critValStepRTab}, i.e. by calling \code{remove(critValStepRTab, envir = envir)}, with \code{envir} the used environment, and from the file system by deleting the corresponding subfolder, i.e. by calling \preformatted{unlink(file.path(R.cache::getCacheRootPath(), dirs), recursive = TRUE),} with \code{dirs} the corresponding subdirectory.
}
}

\seealso{\code{\link{monteCarloSimulation}}, \code{\link{penalty}}, \code{\link{parametricFamily}}, \code{\link{intervalSystem}}, \code{\link{stepFit}}, \code{\link{computeBounds}}}

\references{
Frick, K., Munk, A., Sieling, H. (2014) Multiscale change-point inference. With discussion and rejoinder by the authors. \emph{Journal of the Royal Statistical Society, Series B} \bold{76}(3), 495--580.

Pein, F., Sieling, H., Munk, A. (2017) Heterogeneous change point inference. \emph{Journal of the Royal Statistical Society, Series B}, \bold{79}(4), 1207--1227.
} 

\examples{
# vector of critical values
qVector <- critVal(100L, alpha = 0.5)
# global quantile
qValue <- critVal(100L, alpha = 0.5, output = "value")

# vector can be computed from the global quantile
identical(critVal(100L, q = qValue), qVector)

# for a conservative significance level, stronger confidence statements
critVal(100L, alpha = 0.05)
critVal(100L, alpha = 0.05, output = "value")

# higher significance level for larger detection power, but less confidence
critVal(100L, alpha = 0.99)
critVal(100L, alpha = 0.99, output = "value")
\donttest{
# different parametric family, different intervalSystem, a subset of lengths,
# different penalty and given weights
q <- critVal(100L, alpha = 0.05, family = "hsmuce", intervalSystem = "dyaLen",
             lengths = c(2L, 4L, 16L, 32L), penalty = "weights",
             weights = c(0.4, 0.3, 0.2, 0.1))

# vector of critical values can be given by a vector of length n
vec <- 1:100
vec[c(2L, 4L, 16L, 32L)] <- q
attr(vec, "n") <- 128L
identical(critVal(100L, q = vec, family = "hsmuce", intervalSystem = "dyaLen",
                  lengths = c(2L, 4L, 16L, 32L)), q)

# with a given monte-Carlo simulation for nq = 128 observations
stat <- monteCarloSimulation(128)
critVal(n = 100L, alpha = 0.05, stat = stat)

# the above calls saved and (attempted to) load Monte-Carlo simulations and
# simulated them for nq = 128 observations
# in the following call no saving, no loading and simulation for n = 100
# observations is required, progress of the simulation will be reported
critVal(n = 100L, alpha = 0.05, messages = 1000L,
        options = list(simulation = "vector", load = list(), save = list()))

# only type "vector" will be saved and loaded in the workspace
critVal(n = 100L, alpha = 0.05, messages = 1000L,
        options = list(simulation = "vector", load = list(workspace = "vector"),
                       save = list(workspace = "vector")))

# simulation of type "matrix" will be saved in a RDS file and in a variable
# saving of type "vector" is disabled by passing "",
# different seed is set and number of simulations is reduced to r = 1e3
# to allow faster computation at the price of a less precise result
file <- tempfile(pattern = "file", tmpdir = tempdir(), fileext = ".RDS")
critVal(n = 100L, alpha = 0.05, seed = 1, r = 1e3,
        options = list(simulation = "matrix", load = list(),
          save = list(variable = c("", "test"), RDSfile = c("", file))))
identical(test, monteCarloSimulation(100L, seed = 1, r = 1e3))
identical(readRDS(file), monteCarloSimulation(100L, seed = 1, r = 1e3))
}}

\keyword{nonparametric}
