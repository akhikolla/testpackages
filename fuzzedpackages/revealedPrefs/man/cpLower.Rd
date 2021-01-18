\name{cpLower}
\alias{cpLower}
\alias{print.lowerBound}
\alias{summary.lowerBound}

\title{
Lower bound on the number of GARP-consistent subpopulations.
}
\description{
The \code{cpLower} function computes a Crawford-Pendakur type lower bound on the number of GARP-consistent subpopulations by creating a set of pairwise GARP-violating observations.
}
\usage{
cpLower(x, p, times= 1, afriat.par= 1)
\method{print}{lowerBound}(x, ...)
\method{summary}{lowerBound}(object, ...)
}
\arguments{
  \item{x}{data frame or matrix containing the observed quantities, where each row corresponds to an observation and the columns are types of goods, or an object of class \code{lowerBound} to be used with \code{print},}
  \item{p}{data frame or matrix (of same dimensions as x) containing the corresponding prices,}
  \item{times}{number of times the algorithm is run (the final result is the best of \code{times} results, ie. highest number of clusters found),}
  \item{afriat.par}{the Afriat parameter, a real number in [0,1], which allows a certain level of error in the optimization of choices ; default is 1, ie. no optimization error allowed,}
  \item{object}{object of class \code{lowerBound} as returned by \code{cpLower},}
  \item{...}{additional arguments passed to the \code{print} and \code{summary} methods (unused).}
}
\value{
\code{cpLower} returns an object of class \code{lowerBound} which contains the following elements:
  \item{\code{violators}}{numeric vector containing the indices of observations that are pairwise GARP-inconsistent,}
  \item{\code{n.clust}}{lower bound on the number of types,}
  \item{\code{hist.n.types}}{numeric vector containing the history of numbers of types found during multiple runs of the algorithm.}
  \item{\code{n.types}}{lower bound on the number of types,}
  \item{\code{afriat.par}}{Afriat parameter used in the algorithm.}
}
\details{
For each run of the algorithm, a random permutation of the observations is drawn, and one by one each observation is pairwise-tested against all previously found violators. If the current observation is found pairwise-inconsistent with all previously found violators it is added to 
the set of violators.
}
\references{
Crawford, I. and Pendakur, K. (2013). How many types are there?
\emph{The Economic Journal}, 123(567):77-95.
}
\author{
Julien Boelaert \email{jubo.stats@gmail.com}
}
\seealso{See \link{cpUpper} for the upper bound on the number of types.}
\examples{
# Lower bound for GARP-violating data:
data(noGarp)
cp.low <- cpLower(noGarp$x, noGarp$p)
cp.low
cp.low$violators
}

\keyword{methods}
