\name{prepareData}
\alias{prepareData}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Prepare data and controls for smoothing
}
\description{
This sorts the data, identifies parameters and function value (response), identifies pairs of response values for identical parameter values, and may set some global controls in \code{blackbox.options()}.
}
\usage{
prepareData(data, ParameterNames=NULL, respName=NULL,
            verbose=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
  A data frame including variables in \code{ParameterNames} and \code{respName}  }
  \item{ParameterNames}{
  Names of the variables to be used as predictors of the smoothed surface. If NULL, all columns except the last are assumed to hold parameter values. }
  \item{respName}{
  Name of the variable to be used as response of the smoothed surface. If NULL, the last column is assumed to hold function values.}
 \item{verbose}{
 Whether to print some information (in particular a message if replicate responses values are identical for given parameter values, whichwill be suspect in some applications)}
}
\value{
A data frame with the required variables, ordered by increasing values as in \code{do.call(order,data)}.
This may set some global controls in \code{blackbox.options()} as a side effect.
}
\examples{
require(spaMM)
data(blackcap) ## use dataset as template
sorted_etc <- prepareData(data=blackcap,ParameterNames=c("longitude", "latitude"),
                    respName="means")
}
