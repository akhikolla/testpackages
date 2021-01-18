\name{checkGarp}
\alias{checkGarp}
\alias{checkWarp}
\alias{checkSarp}
\alias{print.axiomTest}
\alias{summary.axiomTest}

\title{
Non-parametric tests of rationality axioms (WARP, SARP, GARP).
}
\description{
Functions for non-parametric tests of compliance with rationality axioms.

Functions \code{checkWarp}, \code{checkSarp} and \code{checkGarp} perform exact tests of the Weak Axiom of Revealed Preferences (WARP), the Strong Axiom of Revealed Preferences (SARP) and the Generalized Axiom of Revealed Preferences (GARP) respectively, to check whether the given data are consistent with utility maximization.
}
\usage{
checkWarp(x, p, afriat.par= 1)
checkSarp(x, p, afriat.par= 1, method= c("deep", "floyd"))
checkGarp(x, p, afriat.par= 1, method= c("floyd"))
\method{print}{axiomTest}(x, ...)
\method{summary}{axiomTest}(object, ...)
}
\arguments{
  \item{x}{data frame or matrix containing the observed quantities, where each row corresponds to an observation and the columns are types of goods, or an object of class \code{axiomTest} for use with \code{print},}
  \item{p}{data frame or matrix (of same dimensions as x) containing the corresponding prices,}
  \item{afriat.par}{the Afriat parameter, real number in [0,1], which allows a certain level of error in the optimization of choices; default is 1, ie. no optimization error allowed,}
  \item{method}{character string: \code{"deep"} for depth-first search (default for Sarp), \code{"floyd"} for Floyd-Warshall algorithm (can be very slow for large datasets when no violations are present),}
  \item{object}{an object of class \code{axiomTest} as returned by eg. \code{checkWarp},}
  \item{...}{additional arguments passed to the \code{print} and \code{summary} methods (unused).}
}
\value{
Functions \code{checkWarp}, \code{checkSarp} and \code{checkGarp} return an object of class \code{axiomTest} which may contain the following elements:
  \item{\code{violation}}{logical value, \code{TRUE} if a violation was found,}
  \item{\code{path}}{last path taken during depth-first search,}
  \item{\code{path.strict}}{(for \code{deep} method) vector of logical values indicating the strictness of direct preference relations in the last path taken,}
  \item{\code{violators}}{vector of indices of axiom-violating observations (only the first violators that were found are reported, this is not a complete list),}
  \item{\code{strict}}{vector of logical values indicating whether the preferences between reported violators (or path) are strict or not,}
  \item{\code{direct.violation}}{logical value, \code{TRUE} if a violation was found in direct preferences,}
  \item{\code{type}}{the violated or non-violated rationality axiom,} 
  \item{\code{method}}{method used for the non-parametric test,}
  \item{\code{afriat.par}}{Afriat parameter used in the algorithm.}
}
\details{
Rationality axioms can be summarized as follows:

WARP: if X directly prefered to Y and X is not equal to Y, then Y cannot be directly prefered to X (WARP is a necessary condition for the existence of a single-valued utility-maximizing demand function consistent with the data).

SARP: if X (in)directly prefered to Y and X is not equal to Y, then Y cannot be (in)directly prefered to X (SARP is a necessary and sufficient condition for the existence of a single-valued utility-maximizing demand function consistent with the data).

GARP: if X (in)directly prefered to Y, then Y cannot be strictly directly prefered to X (GARP is a necessary and sufficient condition for the existence of a (possibly multiple-valued) utility-maximizing demand function consistent with the data). 

If WARP or GARP are violated, then SARP is also violated.
  
Testing of WARP is straightforward by pairwise comparison of revealed preferences. GARP is tested using the Floyd-Warshall algorithm to find the transitive closure of the direct preference relationship. SARP can be tested either using the Floyd-Warshall algorithm, or using a depth-first search that systematically explores the preference relationship in search of a violating cycle.
}
\references{
Varian, H. R. (1982) The Nonparametric Approach to Demand Analysis, 
\emph{Econometrica}, 50(4):945-973.

Varian, H. R. (1984) \emph{Microeconomic Analysis}. New York/London:
Norton, 2nd edition, pp 141-143.
}
\author{
Julien Boelaert \email{jubo.stats@gmail.com}
}

\examples{
# Check axioms for GARP-violating data:
data(noGarp)
summary(checkGarp(noGarp$x, noGarp$p))

# Check axioms for SARP-consistent data:
data(okSarp)
summary(checkSarp(okSarp$x, okSarp$p))
}

\keyword{methods}
