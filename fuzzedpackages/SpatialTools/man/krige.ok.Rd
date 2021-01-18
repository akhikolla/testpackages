\name{krige.ok} 

\alias{krige.ok}

\title{
	Performs Ordinary Kriging
} 

\description{
	Performs Ordinary Kriging using \code{y}, the \eqn{n \times 1} matrix of observed responses,
	\code{V}, the (positive definite) covariance matrix of the 
	observed responses, \code{Vp}, the 
	\eqn{np \times np}
	covariance matrix of the responses to be predicted, and \code{Vop}, 
	the \eqn{n \times np} matrix of covariances between the observed 
	responses and the responses to be predicted.
} 

\usage{
	krige.ok(y, V, Vp, Vop, nsim = 0, Ve.diag = NULL, method = "eigen")
} 

\arguments{
	\item{y}{The vector of observed responses.
		Should be a matrix of size \eqn{n \times 1} or a vector of 
		length \eqn{n}.}
	\item{V}{The covariance matrix of the observed responses.  
		The size is \eqn{n \times n}.} 
	\item{Vp}{The covariance matrix of the responses to be predicted.
		The size is \eqn{np \times np}.} 
	\item{Vop}{The cross-covariance between the observed responses 
		and the responses to be predicted.  The size is 
		\eqn{n \times np}} 
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
	\item{coeff}{A vector of length \eqn{k} containing the estimated
		regression coefficients.} 
	\item{vcov.coeff}{A \eqn{k \times k} matrix containing the
	(estimated) covariance matrix of estimated the regression coefficients.} 
	\item{simulations}{An \eqn{n \times nsim} matrix containing the \code{nsim} realizations of the conditional realizations.  Each column of the matrix represents a realization of the conditional normal distribution.}
If \code{nsim} > 0, this list has class "krigeConditionalSample".
}

\details{
It is assumed that there are \eqn{n} observed data values
and that we wish to make predictions at \eqn{np} locations.

If doing conditional simulation, the Cholesky decomposition should not work when there are coincident locations between the observed data locations and the predicted data locations.    Both the Eigen and Singular Value Decompositions should work.

If user specifies \code{nsim} to be a positive integer, then \code{nsim} conditional realizations of the predictive distribution will be generated.  If this is less than 1, then no conditional simulation is done.  If \code{nsim} is a positive integer, then \code{Ve.diag} must also be supplied.  \code{Ve.diag} is should be a vector of length \eqn{n} specifying the measurement error variances of the observed data.  This information is only used for conditional simulation, so this argument is only needed when \code{nsim} > 0.  When conditional simulation is desired, then the argument \code{method} can be to specify the method used to decompose \code{V}.  Options are "eigen", "chol", or "svd" (Eigen decomposition, Cholesky decomposition, or Singular value decomposition, respectively).  This information is only used for conditional simulation, so this argument is only applicable when \code{nsim} > 0.
}
\references{ Statistical Methods for Spatial Data Analysis, Schabenberger and Gotway (2003).  See p. 226-228.} 

\author{Joshua French} 

\keyword{kriging}
\keyword{ordinary kriging}
\keyword{ordinary} 

\examples{
# create observed and predicted coordinates
ocoords <- matrix(runif(100), ncol = 2)
pcoords <- matrix(runif(200), ncol = 2)

# include some observed locations in the predicted coordinates
acoords <- rbind(ocoords, pcoords)

# create covariance matrix
C3 <- cov.sp(coords = ocoords, sp.type = "matern", sp.par = c(2, 1), smoothness = 1, 
	finescale = 0, error = 0.5, pcoords = acoords)

# generate data with error
y <- rmvnorm(nsim = 1, mu = rep(2, 50), V = C3$V) + rnorm(50, sd = sqrt(.5))

# use universal kriging to make predictions.  Do not do conditional simulation
krige.obj <- krige.ok(as.vector(y), V = C3$V, Vp = C3$Vp, Vop = C3$Vop, 
	nsim = 0)

#Do conditional simulation
krige.obj2 <- krige.ok(as.vector(y), V = C3$V, Vp = C3$Vp, Vop = C3$Vop, 
	nsim = 100, Ve.diag = rep(.5, 50), method = "eigen")

}