\name{summary.cfc.pbasis}
\alias{summary.cfc.pbasis}
\alias{plot.summary.cfc.pbasis}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Summarizing probability-denominated CFC objects
}
\description{
\code{summary} method for class \code{\link{cfc.pbasis}}.
}
\usage{
\method{summary}{cfc.pbasis}(object, ...)
\method{plot}{summary.cfc.pbasis}(x, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{An object of class "cfc.pbasis", usually the result of a call to \code{\link{cfc.pbasis}}.}
  \item{x}{An object of class "summary.cfc.pbasis", usually the result of a call to \code{summary.cfc.pbasis}.}
  \item{...}{Further arguments to be passed to/from other methods.}
}
\value{
The function \code{summary.cfc.pbasis} calculates the average of cumulative incidence and event-free probability functions at each time point across all elements of the object list. If the object is a matrix, it is returned without change.
}

\author{Alireza S. Mahani, Mansour T.A. Sharabiani}

\references{
Mahani A.S. and Sharabiani M.T.A. (2019). Bayesian, and Non-Bayesian, Cause-Specific Competing-Risk Analysis for Parametric and Nonparametric Survival Functions: The R Package CFC. Journal of Statistical Software, 89(9), 1-29. doi:10.18637/jss.v089.i09
}

\seealso{
The model fitting function is \code{\link{cfc.pbasis}}. See \code{\link{summary}} and \code{\link{plot}} for descriptions of the generic methods. See \code{\link{cfc.tbasis}} for time-denominated CFC, as well as usage examples.

}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.

