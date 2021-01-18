\name{ojaMedianControl}
\alias{ojaMedianControl}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Tuning Parameters for the Function 'ojaMedian'}
\description{
 Tuning parameters for the algorithms used by function \code{\link{ojaMedian}.}
}

\usage{
ojaMedianControl(sigmaInit = 0, sigmaAda = 20, adaFactor = 0.5, 
                  iter = 1e+06, useAllSubsets = FALSE, 
                  nSubsetsUsed = 1000, sigmaLog10Dec = 10, 
                  storeSubDet = TRUE, eps = 0.1, chi2 = 0.95, 
                  samples = 20, maxlines = 1000, S1 = cov, 
                  S2 = cov4, S1args = list(), S2args = list(),
                  volume = 1e-6, boundedExact = T, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{sigmaInit}{(for the \code{evo} algorithm): Set the initial variance of the mutation vector in the first run.}
  \item{sigmaAda}{(for the \code{evo} algorithm): Defines after how many mutations the variance of the mutation vector is adjusted. }
  \item{adaFactor}{(for the \code{evo} algorithm): Defines the level of adjustment of the mutation vector.}
  \item{iter}{(for the \code{evo} algorithm): The maximum number of iterations. If the algorithm does not converge, it stops after \code{iter} - iterations.}
  \item{useAllSubsets}{(for the \code{evo} algorithm): A logical flag. If it is set all datapoints and resulting simplices are taken into account for the calculation.}
  \item{nSubsetsUsed}{(for the \code{evo} algorithm): If useAllSubsets is not set, this determines how many, randomly selected, datapoints are taken into account. }
  \item{sigmaLog10Dec}{(for the \code{evo} algorithm): This is an abort criterion. If the logarithmised initial variance differs more than sigmaLog10Dec from the actual, logarithmised variance, the algorithm stops.}
  \item{storeSubDet}{(for the \code{evo} algorithm): A boolean flag. If it is set subdeterminants are stored. This should always been set to \code{TRUE} if \eqn{6*(dim -1)*nSubsetsUsed<CPU-Cache}. }
  \item{eps}{(for the \code{grid} algorithm): This is the abort criterion. If the grid becomes denser than this threshold the algorithm stops.}
  \item{chi2}{(for the \code{grid} algorithm): This is the test niveau of the test, if a grid point could be used as a Oja-Median or not.}
  \item{samples}{(for the \code{grid} algorithm): This determines how many additional hyperplanes are taken after every run.}
  \item{maxlines}{(for the \code{exact} algorithm): This determines how many intersection lines are investigated in addtion to the one with the steepest gradient.}
  \item{S1}{(for the \code{evo} and \code{grid} algorithms): Passed on to \code{\link[ICS]{ics}} to compute the invariant coordinate system. Default is \code{\link{cov}}.}
  \item{S2}{(for the \code{evo} and \code{grid} algorithms): Passed on to \code{\link[ICS]{ics}} to compute the invariant coordinate system. Default is \code{\link[ICS]{cov4}}.}
  \item{S1args}{(for the \code{evo} and \code{grid} algorithms): Optional arguments for \code{S1} passed on to \code{\link[ICS]{ics}} to compute the invariant coordinate system.}
  \item{S2args}{(for the \code{evo} and \code{grid} algorithms): Optional arguments for \code{S2} passed on to \code{\link[ICS]{ics}} to compute the invariant coordinate system.}
  \item{volume}{(for the \code{bounded_exact} algorithm):  is the desired size of the bounded region, which is selected as a part of the original volume. Here the volume is calculated as the volume of a minimal multivariate circumscribed rectangle with edges parallel to the coordinate axes.}
  \item{boundedExact}{(for the \code{bounded_exact} algorithm):  setting this parameter to \code{FALSE} stops the algorithm after the bounded region is found, and its center is reported as an approximation of the median.}
  \item{...}{Additional hidden option}
}

\author{Daniel Fischer and Klaus Nordhausen}

\seealso{\code{\link{ojaMedian}}, also for references and examples.  }

\examples{
## Show the default settings:
str(ojaMedianControl())
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{misc}

