\name{summary.cfc.survreg}
\alias{summary.cfc.survreg}
\alias{plot.summary.cfc.survreg}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Summarizing and plotting output of \code{cfc.survreg}
}
\description{
\code{summary} and method for class \code{\link{cfc.survreg}}.
}
\usage{
\method{summary}{cfc.survreg}(object, obs.idx = "all", ...)
\method{plot}{summary.cfc.survreg}(x, which = c(1, 2), ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{An object of class "cfc.survreg", usually the result of a call to \code{\link{cfc.survreg}}.}
  \item{obs.idx}{Index of observations to calculate mean cumulative incidence for; defaults to all observation.}
  \item{...}{Further arguments to be passed to/from other methods.}
  \item{x}{An object of class \code{summary.cfc.survreg}, usually the output of \code{summary.cfc.survreg}.}
  \item{which}{Vector of integers, indicating which plot(s) must be produced: 1) cumulative incidence functions, one per cause, as a function of time-to-index, all in the same plot, 2) comparison of cumulative incidence function with/without competing-risk adjustment. The unadjusted figure is equivalent to 1 minus the Kaplan-Meyer (i.e., survival) function.}
}
\value{
\code{summary.cfc.surveeg} produces a matrix of dimensions \code{length(object$tout)} (number of time points) by \code{object$K} (number of causes). See description of \code{which} aregument for \code{plot.summary.cfc.survreg}.
}

\author{Mansour T.A. Sharabiani, Alireza S. Mahani}

\references{
Mahani A.S. and Sharabiani M.T.A. (2019). Bayesian, and Non-Bayesian, Cause-Specific Competing-Risk Analysis for Parametric and Nonparametric Survival Functions: The R Package CFC. Journal of Statistical Software, 89(9), 1-29. doi:10.18637/jss.v089.i09
}

\seealso{
The model fitting function is \code{\link{cfc.survreg}}. See \code{\link{summary}} and \code{\link{plot}} for descriptions of the generic methods. For more flexible ways of cause-specific competing-risk analysis, see \code{\link{cfc}}.

}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.

