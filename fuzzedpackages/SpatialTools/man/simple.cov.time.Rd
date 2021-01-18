\name{simple.cov.time} 

\alias{simple.cov.time}

\title{
	Calculates temporal covariance based on distance matrix
} 

\description{
	Calculates a temporal covariance using a (Euclidean) distance matrix \code{T}.  Not intended to be used directly by user (though it may be helpful to some).  It is used in the \code{covets} function.  No argument or error checking is provided for this function.
} 

\usage{
	simple.cov.time(T, t.type, t.par)
} 

\arguments{
	\item{T}{A distance matrix.} 
	\item{t.type}{A character vector specifying the temporal covariance type.  Only "ar1" is currently implemented.} 
	\item{t.par}{A vector of length 1 specifying the strength of dependence of the covariance function.} 
} 

%\details{
%}

\value{
	Returns a covariance matrix. 
}


%\references{ R core development. "Writing R Extentions". 2007. } 

\author{Joshua French} 

\keyword{covariance} 
\keyword{spatial} 

\seealso{ ~ \code{cov.st} } 

\examples{
	T <- dist1(matrix(1:10))
	simple.cov.time(T = T, t.type = "ar1", t.par = .5)
}