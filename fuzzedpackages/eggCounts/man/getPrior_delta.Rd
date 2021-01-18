\name{getPrior_delta}
\alias{getPrior_delta}
\title{Get prior parameters from Beta distribution}
\usage{
getPrior_delta(lower, upper, p = 0.7, mode, conc, plot = TRUE)
}
\arguments{
\item{lower,upper,p}{numeric. Prior belief about the reduction. There is \code{p} probability that the reduction is between \code{lower} and \code{upper}. Not used if \code{mode} and \code{conc} are supplied.}

\item{mode,conc}{numeric. Prior belief about the reduction. The mode and concentration parameters of a beta distribution. Higher concentration indicates smaller variance. Not used if \code{lower} and \code{upper} thresholds are supplied.}

\item{plot}{logical. If TRUE, the prior distribution is plotted after parameters are found.}
}
\value{
Returns Beta prior parameters for \eqn{\delta} and the printed argument to use in a \code{\link{fecr_stan}()} or a \code{\link{fec_stan}()} function call.
}
\description{
Compute the shape parameters from a Beta distribution for \eqn{\delta} based on some prior belief.
}
\author{
Tea Isler \cr
Craig Wang 
}
\details{
The \code{\link[rootSolve]{multiroot}} function from \strong{rootSolve} package is used to compute the parameters.
}
\examples{
# there is 80\% probability that the reduction is between 60\% and 90\%
getPrior_delta(lower = 0.6, upper = 0.9, p = 0.8)
}
