\name{cov.sp} 

\alias{cov.sp}

\title{
	Calculates spatial covariance
} 

\description{
	Calculates spatial covariance matrix of the observed responses, and 
	possibly, the responses to be predicted.  
	If \code{pcoords} is not provided, then only \code{V}, 
	the covariance matrix of the observed responses will be returned.  
	If \code{pcoords} is provided, then \code{Vp} and \code{Vop} 
	(the covariance matrix for predicted responses and between observed and 
	predicted responses, respectively) will also be returned.	
} 

\usage{
	cov.sp(coords, sp.type = "exponential", 
		sp.par = stop("specify sp.par argument"), 
		error.var = 0, smoothness = 0.5, finescale.var = 0, 
		pcoords = NULL, D = NULL, Dp = NULL, Dop = NULL)
} 

\arguments{
	\item{coords}{A numeric matrix of size \eqn{n \times d} containing the observed data locations.} 
	\item{sp.type}{A character vector specifying the spatial covariance type.  Valid types are currently exponential, gaussian, matern, matern2, and spherical.} 
	\item{sp.par}{A vector of length 2 specifying the scale and strength of dependence of the covariance function.  The first element is the variance of the underlying spatial process (also known as the hidden or latent spatial process).  This value is also called the partial sill.   The second element is the strength of dependence between responses.}
	\item{error.var}{A non-negative number indicating the variance of the error term.} 
	\item{smoothness}{A positive number indicating the variance of the error term.} 
	\item{finescale.var}{A non-negative positive number indicating the finescale variability.  The is also called the microscale variance} 
	\item{pcoords}{A numeric matrix of size \eqn{np \times d} containing the locations of the responses to be predicted.} 
	\item{D}{The Euclidean distance matrix for the \code{coords} matrix.  Must be of size \eqn{n \times n}.}
	\item{Dp}{The Eucliean distance matrix for the \code{pcoords} matrix.  Must be of size \eqn{np \times np}.}
	\item{Dop}{The Euclidean intersite distance matrix between the locations in \code{coords} and the locations in \code{pcoords}.  Must be of size \eqn{n \times np}.}
} 

\details{
  The spatial covariance functions are parameterized in a manner 
  consistent with the \code{cov.spatial} function in the 
  \code{geoR} package.  The \code{matern2} covariance function is an alternative covariance function suggested by Handcock and Wallis (1994).  
  The benefit of this parameterization is that the range parameter is that it allows the effective range to be less dependent on the smoothness parameter.

	The \code{D}, \code{Dp}, and \code{Dop} arguments are supplied to decrease the number of necessary computations needed when performing repetitive analysis or simulations.  It is probably in the user's interest to not supply these arguments unless the duration of analysis is an important consideration.  Note that these arguments override the information given in \code{coords} and \code{pcoords}, i.e., if dist1(coords) != D, then D is used in subsequent calculations, etc.  This could create problems.
}

\value{
	Returns a list with the following elements:
\item{V}{ The covariance matrix for the observed responses.   Will be of size \eqn{n \times n}.}
\item{Vp}{ The covariance matrix for the predicted responses.  Only returned if \code{pcoords} is supplied.  Will be of size \eqn{np \times np}.}
\item{Vp}{ The covariance matrix between the observed responses and the predicted responses.  Only returned if \code{pcoords} is supplied.  Will be of size \eqn{n \times np}.}
}


\references{ M.S. Handcock, J.R. Wallis.  An approach to statistical spatial-temporal modeling of meteorological fields (with discussion).  
Journal of the American Statistical Association, 89 (1994), pp. 368--390. } 

\author{Joshua French} 

\keyword{covariance} 
\keyword{spatial} 

\seealso{\code{simple.cov.sp} } 

\examples{
    coords <- matrix(rnorm(30), ncol = 3)
    cov.sp(coords = coords, sp.type = "exponential", sp.par = c(2, 1),
        error.var = 1)
}