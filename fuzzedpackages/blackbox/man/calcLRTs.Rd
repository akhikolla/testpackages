\name{calcLRTs}
\alias{calcLRTs}
\title{
Compute (profile) likelihood ratio tests
}
\description{
Assuming that \code{\link{calcPredictorOK}} and \code{\link{maximizeOK}} have been first run and that the predicted response surface is a likelihood surface , this performs likelihood ratio (LR) tests for a list of parameter points. Profiles are computed if appropriate, i.e. is the point is lower-dimensional than the the parameter space.
}
\usage{
calcLRTs(testPointList, cleanResu = "")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{testPointList}{
A list of points in predictor (parameter) space. Each point is a numeric vector or list with named elements , the names being those of some parameters.}
  \item{cleanResu}{
  A connection, or a character string naming a file for some nicely formated output. If \code{""} (the default), print to the standard output connection.
}
}
\value{
Return a list with information about each LR test, except for tests that could not be performed (e.g. if the tested point is ousdie of the convex envelope of the parameter points from which the predictor has been built). The names of this list's elements are constructed from the tested points.
Eachelement is itself a list with elements
\item{LRT}{The LR statistics (twice the differnece in log-likelihood between maximized likelihood and profile for the input parameters)}
\item{pval}{Associated Pvalue by standard chi-square approximation}
\item{profpt}{Information about the profile point for the input parameters}
\item{maxpt}{Information about the maximum likelihood point}
and other elements, not documented here.
}
