\name{createModel}
\alias{createModel}
\title{Create and Initialize the Rcpp_bei Model Object}
\description{
createModel creates and initializes an \code{\link{Rcpp_bei}} object.
}
\usage{
createModel(data.obj,K)
}
\arguments{
	\item{data.obj}{\code{\link{EditIn.data}} object generated from \code{\link{readData}}.}
	\item{K}{maximum number of DP mixture components.}
}
\value{
createModel returns an \code{\link{Rcpp_bei}} model-object. The returned model object will be referenced in all subsequent calls.
}
\seealso{\code{\link{Rcpp_bei}}} 