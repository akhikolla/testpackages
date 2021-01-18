\name{summary.cfc.tbasis}
\alias{summary.cfc.tbasis}
\alias{plot.summary.cfc.tbasis}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Summarizing time-denominated CFC objects
}
\description{
\code{summary} method for class \code{\link{cfc.tbasis}}.
}
\usage{
\method{summary}{cfc.tbasis}(object,
  MARGIN = if (class(object)[2] == "matrix") NULL else 1,
  ...)
\method{plot}{summary.cfc.tbasis}(x, t = 1, ci = 0.95, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{An object of class "cfc.tbasis", usually the result of a call to \code{\link{cfc.tbasis}}.}
  \item{MARGIN}{Dimensions of cumulative incidence and event-free probability arrays to keep while averaging the remaining dimensions. If the \code{cfc.tbasis} object is a matrix, no averaging is performed, and the function returns the object unchanged. Note that for list objects, \code{MARGIN} must include \code{1} since it is meaningless to average out the time/probability dimension (which is always the first one).}
  \item{x}{An object of class "summary.cfc.tbasis", usually the result of a call to \code{summary.cfc.tbasis}.}
  \item{t}{Regular time step, or vector of time values, used for producing cumulative incidence and event-free probability plots.}
  \item{ci}{Confidence interval used in cumulative incidence and event-free probability plots.}
  \item{...}{Further arguments to be passed to/from other methods.}
}
\value{
The function \code{summary.cfc.tbasis} calculates the average of cumulative incidence and event-free probability functions as directed by \code{MARGIN}. For example, if the element \code{ci1} of the object list is three-dimensional, then using \code{MARGIN=c(1,2)} causes the last dimension to be averaged out.
}

\author{Alireza S. Mahani, Mansour T.A. Sharabiani}

\references{
Mahani A.S. and Sharabiani M.T.A. (2019). Bayesian, and Non-Bayesian, Cause-Specific Competing-Risk Analysis for Parametric and Nonparametric Survival Functions: The R Package CFC. Journal of Statistical Software, 89(9), 1-29. doi:10.18637/jss.v089.i09
}

\seealso{
The model fitting function is \code{\link{cfc.tbasis}}. See \code{\link{summary}} and \code{\link{plot}} for descriptions of the generic methods. See \code{\link{cfc.pbasis}} for probability-denominated CFC, as well as usage examples.

}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.

