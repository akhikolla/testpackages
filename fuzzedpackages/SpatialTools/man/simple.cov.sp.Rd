\name{simple.cov.sp} 

\alias{simple.cov.sp}

\title{
	Calculates spatial covariance based on distance matrix
} 

\description{
	Calculates a spatial covariance using a (Euclidean) distance matrix \code{D}.  Not intended to be used directly by user (though it may be helpful to some).  It is strongly recommended that you use the \code{cov.sp} function.  No argument or error checking is provided for this function.
} 

\usage{
	simple.cov.sp(D, sp.type, sp.par, error.var, smoothness, finescale.var)
} 

\arguments{
	\item{D}{A distance matrix between locations} 
	\item{sp.type}{A character vector specifying the spatial covariance type.  Valid types are currently exponential, gaussian, matern, and spherical.} 
	\item{sp.par}{A vector of length 2 specifying the scale and dependence of the covariance function.  The first element refers to the variance of the hidden process (sometimes this is called the partial sill) while the second elements determines the strength of dependence between locations.} 
	\item{error.var}{A non-negative number indicating the variance of the error term.} 
	\item{smoothness}{A positive number indicating the variance of the error term.} 
	\item{finescale.var}{A non-negative positive number indicating the finescale variability.  The is also called the microscale variance} 
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

\seealso{ ~ \code{cov.sp} } 

\examples{
	coords <- matrix(rnorm(30), ncol = 3)
	D <- dist1(coords)
	simple.cov.sp(D = D, sp.type = "exponential", sp.par = c(2, 1),
		error.var = 1, smoothness = 0.5, finescale.var = 0)
}