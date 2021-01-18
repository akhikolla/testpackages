\name{cfc.prepdata}
\alias{cfc.prepdata}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Utility function for CFC data preparation
}
\description{
Preparing a data frame and formulas for cause-specific competing-risk survival analysis. It expands the multi-state status column into a series of binary columns by treating an event for a cause as censoring for all other causes.
}
\usage{
cfc.prepdata(formul, dat)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formul}{Original survival formula.}
  \item{dat}{Original data frame, with \code{status} column being an integer with values from \code{0} to \code{K}. The value \code{0} represents right-censoring, while \code{1} to \code{K} represent the \code{K} mutually-exclusive events.}
}
\details{
The output data frame will have \code{K} new binary status columns. The \code{K} new status columns will be named \code{"status_1"}, \code{"status_2"} through \code{"status_<K>"}. Each of the output formulas in \code{formula.list} field will have the corresponding \code{status}. Column "status_1" will be \code{1} wherever \code{status} equals \code{1} in original data frame, and \code{0} elsewhere, and similarly for the remaining \code{K-1} newly-added status columns.
}
\value{
A list with the following elements:
\item{K}{Number of causes.}
\item{dat}{Expanded data frame.}
\item{formula.list}{A list of \code{K} formulas, each corresponding to one of the cause-specific survival models to be estimated. See details.}
\item{formula.noresp}{A formula with no left-hand side (time and status variables). This can be used for preparing the model matrix for prediction data sets, which can possibly have no response.}
\item{tmax}{Maximum time to event/censoring extracted from original data frame. This can be used, e.g., during competing-risk analysis.}
}
\author{
Mansour T.A. Sharabiani, Alireza S. Mahani
}

\references{
Mahani A.S. and Sharabiani M.T.A. (2019). Bayesian, and Non-Bayesian, Cause-Specific Competing-Risk Analysis for Parametric and Nonparametric Survival Functions: The R Package CFC. Journal of Statistical Software, 89(9), 1-29. doi:10.18637/jss.v089.i09
}

%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{
data(bmt)
prep.out <- cfc.prepdata(Surv(time, cause) ~ platelet + age + tcell, bmt)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
