\name{calcGCV}
\alias{calcGCV}
\title{
Estimate smoothing parameters by generalized cross-validation (GCV)
}
\description{
  Smoothing is based on prediction in a linear mixed model (\dQuote{Kriging}) with non-zero residual variance. The correlation function for the random effect is the Matern function with argument the Euclidian distance between scaled coordinates (x/scale). The Matern function also has a smoothness parameter. These parameters are by default estimated by GCV. For large data sets (say >2000 rows), it is strongly recommended to select a subset of the data using \code{GCVptnbr}, as GCV will otherwise be very slow.

}
\usage{
calcGCV(sorted_data=data, data, CovFnParam = NULL, GCVptnbr = Inf,
       topmode = FALSE, verbose = FALSE, cleanResu = "",
       force=FALSE, decreasing=FALSE,
       verbosity = blackbox.getOption("verbosity"),
       optimizers = blackbox.getOption("optimizers"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{sorted_data}{
  A data frame with both predictor and response variance, sorted and with attributes, as produced by \code{\link{prepareData}}
}
  \item{data}{Obsolete, for Migraine back-compatibility, should not be used.}
  \item{CovFnParam}{
  Optional fixed values of scale factors for each predictor variable. Smoothness should not be included in this argument.
}
  \item{GCVptnbr}{
  Maximum number of rows selected for GCV.
}
  \item{topmode}{
  Controls the way rows are selected. For development purposes, should not be modified
}
  \item{verbose}{
  Whether to print some messages or not. Distinct from \code{verbosity}
  }
  \item{verbosity}{Distinct from \code{verbose}. See \code{verbosity} in \code{\link{blackbox.options}}}
  \item{cleanResu}{
  A connection, or a character string naming a file for some nicely formated output. If \code{""} (the default), print to the standard output connection.
}
  \item{force}{ Boolean. Forces the analysis of data without pairs of response values for given parameter values. }
  \item{optimizers}{A vector of) character strings, from which the optimization method is selected. Default is \code{nloptr} with its own \code{"NLOPT_LN_BOBYQA"} method. See the source of the function for other methods (the latter being subject to change with little notice).}
  \item{decreasing}{ Boolean. Use TRUE if you want the result to be used in function maximization rather than minimization.}
}
\value{
A list with the following elements
\item{CovFnParam}{Scale parameters \bold{and} smoothness parameter of the Matern correlation function}
\item{lambdaEst}{Ratio of residual variance over random effect variance}
\item{pureRMSE}{Estimate of root residual variance}
and possibly other elements.

Global options \code{CovFnParam} is modified as a side effect.
}
\references{
  Golub, G. H., Heath, M. and Wahba, G. (1979) Generalized Cross-Validation as a method for choosing a good ridge parameter.
  Technometrics 21: 215-223.
}
\examples{
# see example on main doc page (?blackbox)
}
