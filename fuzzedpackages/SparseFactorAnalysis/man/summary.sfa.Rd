\name{summary.sfa}
\alias{summary.sfa}
\title{
Summaries from sparse factor analysis.
}
\description{
	Returns a summary for the dimensionality estimate, first and second dimensions, and additional dimensions as desired.
	
}
\usage{
\method{summary}{sfa}(object,... )
}
\arguments{
  \item{object}{
Object of type \env{sfa}.
}
\item{...}{Additional items to pass to \env{summary}.  Options below.}

}



\details{
Generates a table for an object of class \env{sfa}.  Additional arguments to pass \env{summary} below.

\env{interval}{
	A number between 0 and 1. Length of symmetric posterior credible interval.}
	
\env{topbottom}{ A positive integer.  How many observations at the top and bottom of each dimension to print.}

\env{print.dims}{ Optional. Additional dimensions besides the first two for which to return the top and bottom \env{topbottom} observations.}
}
\references{
Ratkovic, Marc and Tingley, Dustin.  2015.  "Sparse Estimation with Uncertainty: Subgroup Analysis in Large Dimensional Design."  Working paper.
}
\seealso{\link{sfa}}
\examples{

##  See sfa() for example.

}
