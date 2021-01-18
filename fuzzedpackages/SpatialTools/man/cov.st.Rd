\name{cov.st} 

\alias{cov.st}

\title{
	Calculates spatio-temporal covariance
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
	cov.st(coords, time, sp.type = "exponential", 
	    sp.par = stop("specify sp.par argument"),
 	    error.var = 0, smoothness = 0.5, finescale.var = 0, 
	    t.type = "ar1", t.par = .5, 
	    pcoords = NULL, ptime = NULL, 
	    D = NULL, Dp = NULL, Dop = NULL, 
	    T = NULL, Tp = NULL, Top = NULL)
} 

\arguments{
	\item{coords}{A numeric matrix of size \eqn{n \times d} containing the observed data locations.} 
	\item{time}{A numeric matrix of size \eqn{n \times 1} containing the times at which the data was observed.} 
	\item{sp.type}{A character vector specifying the spatial covariance type.  Valid types are currently exponential, gaussian, matern, and spherical.} 
	\item{sp.par}{A vector of length 2 specifying the scale and dependence of the covariance function.  The first element refers to the variance of the hidden process (sometimes this is called the partial sill) while the second elements determines the strength of dependence between locations.}
	\item{error.var}{A non-negative number indicating the variance of the error term.} 
	\item{smoothness}{A positive number indicating the variance of the error term.} 
	\item{finescale.var}{A non-negative positive number indicating the finescale variability.  The is also called the microscale variance}
	\item{t.type}{A character vector indicating the temporal dependance structure.  Currently, only "ar1" is implemented.}
	\item{t.par}{A numeric vector of length 1 indicating the strength of temporal dependence.} 
	\item{pcoords}{A numeric matrix of size \eqn{np \times d} containing the locations of the responses to be predicted} 
	\item{ptime}{A numeric matrix of size \eqn{np \times 1} containing the times at which the responses are to be predicted.} 
	\item{D}{The Euclidean distance matrix for the coords matrix.  Must be of size \eqn{n \times n}.}
	\item{Dp}{The Eucliean distance matrix for the pcoords matrix.  Must be of size \eqn{np \times np}.}
	\item{Dop}{The Eucliean intersite distance matrix between the locations in coords and the locations in pcoords.  Must be of size \eqn{n \times np}.}
	\item{T}{The Euclidean distance matrix for the time matrix.  Must be of size \eqn{n \times n}.}
	\item{Tp}{The Eucliean distance matrix for the ptime matrix.  Must be of size \eqn{np \times np}.}
	\item{Top}{The Eucliean intertime distance matrix between the times in \code{time} and \code{ptime}.  Must be of size \eqn{n \times np}.}

} 

\details{
	At this point, this function only implements a separable spatio-temporal covariance funcation.  If \eqn{h} is the distance between
	two sites, and \eqn{t} is the temporal lag between the times when the associated responses were observed, then the covariance function
	\eqn{C(h,t) = Cs(h) \times Ct(t)} where \eqn{Cs} is a spatial covariance function corresponding to the 
	exponential, matern, gaussian, or spherical and \eqn{Ct} is the temporal covariance function corresponding to an ar1 process with 
	\eqn{Ct(t) = \phi ^ t}.

	The \code{D}, \code{Dp}, \code{Dop}, \code{T}, \code{Tp}, \code{Top} arguments are supplied to decrease the number of necessary computations needed when performing repetitive analysis or simulations.  It is probably in the user's interest to not supply these arguments unless the duration of analysis is an important consideration.  Note that these arguments override the information given in \code{coords}, \code{pcoords}, \code{time}, and \code{prime}, i.e., if dist1(coords) != D, then D is used in subsequent calculations, etc.  This could create problems.
}

\value{
	Returns a list with the following elements:
\item{V}{ The covariance matrix for the observed responses. }
\item{Vp}{ The covariance matrix for the predicted responses.  Only returned if \code{pcoords} is supplied.}
\item{Vp}{ The covariance matrix between the observed responses and the predicted responses.  Only returned if \code{pcoords} is supplied.  Will be of size \eqn{n \times np}}
}


%\references{ R core development. "Writing R Extentions". 2007. } 

\author{Joshua French} 

\keyword{covariance} 
\keyword{spatial} 
\keyword{time}
\keyword{spatio-temporal}

\seealso{\code{simple.cov.sp} } 

\examples{
	coords <- matrix(rnorm(30), ncol = 3)
	pcoords <- matrix(rnorm(90), ncol = 3)
	time <- 1:10
	ptime <- 1:30
	cov.st(coords = coords, time = time, sp.type = "exponential", 
	    sp.par = c(2, 1), error.var = 1, t.type = "ar1", t.par = .5, 
	    pcoords = pcoords, ptime = ptime)
}