\name{krige.sk} 

\alias{krige.sk}

\title{
	Performs simple kriging
} 

\description{
	Performs simple kriging using \code{y}, a vector of length \eqn{n},
	\code{V}, the (positive definite) covariance matrix of the 
	observed responses, \code{Vp}, the 
	\eqn{np \times np}
	covariance matrix of the responses to be predicted, \code{Vop}, 
	the \eqn{n \times np} matrix of covariances between the observed 
	responses and the responses to be predicted, and \code{m}, a numeric vector 
  of length 1 identifying the value of the mean
  for each response.
} 

\usage{
	krige.sk(y, V, Vp, Vop, m = 0, nsim = 0, Ve.diag = NULL, method = "eigen")
} 

\arguments{
	\item{y}{The vector of observed responses.
		Should be a matrix of size \eqn{n \times 1} or a vector of 
		length \eqn{n}.}
	\item{V}{The covariance matrix of the observed responses.  
		The size is \eqn{n \times n}.} 
	\item{Vp}{The covariance matrix of the responses to be predicted.
		The size is \eqn{np \times np}} 
	\item{Vop}{The cross-covariance between the observed responses 
		and the responses to be predicted.  The size is 
		\eqn{n \times np}.} 
  	\item{m}{A numeric vector of length 1 giving the mean of each response.}
	\item{nsim}{The number of simulated data sets to sample from the conditional predictive distribution.}
	\item{Ve.diag}{A vector of length \eqn{n} specifying the measure error variances of the observed data.  Only needed if \code{nsim} > 0.}
	\item{method}{The method for decomposing \code{V} in conditional simulation.  Default is "eigen", for the Eigen decomposition.  Alternatives are "chol" (Cholesky) and "svd" (Singular Value Decomposition).}
} 

\value{
The function returns a list containing the following objects:
	\item{pred}{A vector of length \eqn{np} containing the predicted
		responses.}
	\item{mspe}{A vector of length \eqn{np} containing the 
		mean-square prediction error of the predicted
		responses.}
	\item{simulations}{An \eqn{n \times nsim} matrix containing the \code{nsim} realizations of the conditional realizations.  Each column of the matrix represents a realization of the conditional normal distribution.}
	\item{mean}{The mean value (m) originally provided to the function}.  
If \code{nsim} > 0, this list has class "krigeConditionalSample".
}

\details{
It is assumed that there are \eqn{n} observed data values
and that we wish to make predictions at \eqn{np} locations.
The mean is subtracted from each value of \code{y} before determining the kriging weights, 
and then the mean is added onto the predicted response.

If doing conditional simulation, the Cholesky decomposition should not work when there are coincident locations between the observed data locations and the predicted data locations.    Both the Eigen and Singular Value Decompositions should work.

If user specifies \code{nsim} to be a positive integer, then \code{nsim} conditional realizations of the predictive distribution will be generated.  If this is less than 1, then no conditional simulation is done.  If \code{nsim} is a positive integer, then \code{Ve.diag} must also be supplied.  \code{Ve.diag} is should be a vector of length \eqn{n} specifying the measurement error variances of the observed data.  This information is only used for conditional simulation, so this argument is only needed when \code{nsim} > 0.  When conditional simulation is desired, then the argument \code{method} can be to specify the method used to decompose \code{V}.  Options are "eigen", "chol", or "svd" (Eigen decomposition, Cholesky decomposition, or Singular value decomposition, respectively).  This information is only used for conditional simulation, so this argument is only applicable when \code{nsim} > 0.
}
\references{ Statistical Methods for Spatial Data Analysis, Schabenberger and Gotway (2003).  See p. 226-228.} 

\author{Joshua French} 

\keyword{kriging}
\keyword{simple kriging}
\keyword{simple} 

\examples{
	data(toydata)
	y <- as.vector(toydata$y)
	V <- toydata$V
	Vp <- toydata$Vp
	Vop <- toydata$Vop
	krige.sk(y, V, Vp, Vop, m = 2)
}