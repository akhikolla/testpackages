\name{maxlik.cov.sp} 

\alias{maxlik.cov.sp}

\title{
	Determines maximum likelihood estimates of covariance parameters
} 

\description{
	Estimates covariance parameters of spatial covariance functions using maximum likelihood or restricted maximum likelihood.  See \code{cov.sp} for more details of covariance functions to be estimated.  
} 

\usage{
	maxlik.cov.sp(X, y, coords, sp.type = "exponential", 
	    range.par = stop("specify range.par argument"), 
	    error.ratio = stop("specify error.ratio argument"), 
	    smoothness = 0.5,  
	    D = NULL, reml = TRUE, lower = NULL, upper = NULL, 
	    control = list(trace = TRUE), optimizer="nlminb")
} 

\arguments{
	\item{X}{A numeric matrix of size \eqn{n \times k} containing the design matrix of the data locations.} 
	\item{y}{A vector of length \eqn{n} containing the observed responses.} 
	\item{coords}{A numeric matrix of size \eqn{n \times d} containing the locations of the observed responses.} 
	\item{sp.type}{A character vector specifying the spatial covariance type.  Valid types are currently exponential, gaussian, matern, and spherical.} 
	\item{range.par}{An initial guess for the spatial dependence parameter.}
	\item{error.ratio}{A value non-negative value indicating the ratio \code{error.var/sp.par[1]}.} 
	\item{smoothness}{A positive number indicating the smoothness of the matern covariance function, if applicable.} 
	\item{D}{The Euclidean distance matrix for the coords matrix.  Must be of size \eqn{n \times n}.}
	\item{reml}{A boolean value indicating whether restricted maximum likelihood estimation should be used.  Defaults to TRUE.}
	\item{lower}{A vector giving lower bounds for the covariance parameters \code{sp.par[2]}, \code{error.ratio}, and \code{smoothness} (when the model is matern).  Order matters!  If not given defaults to an upper bound of Inf for \code{sp.par[2]}, 1 for \code{error.ratio}, and 10 for \code{smoothness}.}
	\item{upper}{A vector giving upper bounds for the covariance parameters \code{sp.par[2]}, \code{error.ratio}, and \code{smoothness} (when the model is matern).  Order matters!  If not given defaults to an upper bound of Inf for \code{sp.par[2]}, 1 for \code{error.ratio}, and 10 for \code{smoothness}.}
	\item{control}{A list giving tuning parameters for the \code{nlminb} function.  See \code{nlminb} for more details.}
	\item{optimizer}{A vector describing the optimization function to use for the optimization.  Currently, only \code{nlminb} is an acceptable value.}
} 

\value{
	Returns a list with the following elements:
\item{sp.type}{ The covariance form used. }
\item{sp.par}{ A vector containing the estimated variance of the hidden process and the spatial dependence.}
\item{error.var}{ The estimated error variance.}
\item{smoothness}{ The smoothness of the matern covariance function.}
\item{par}{ The final values of the optimization parameters.  Note that these will not necessarily match \code{sp.par}, \code{error.var}, and \code{smoothness} because of the reparameterization.}
\item{convergence}{ Convergence message from \code{nlminb}.}
\item{message}{ Message from \code{nlminb}.}
\item{iterations}{ Number of iterations for optimization to converge.}
\item{evaluations}{ Evaluations from \code{nlminb}.}
}

\details{
	When doing the numerical optimizaiton, the covariance function is reparameterized slightly to speedup computation.  
	Specifically, the variance parameter for the process of interest,\code{sp.par[1]}, is profiled out,
	and the \code{error.var} parameter is parameterized as \code{sp.par[1] * error.ratio}, where \code{error.ratio = error.var/sp.par[1]}.
}

%\references{ R core development. "Writing R Extentions". 2007. } 

\author{Joshua French} 

\keyword{covariance} 
\keyword{spatial} 
\keyword{maximum likelihood}

\seealso{ \code{cov.st} } 

\examples{
#generate 20 random (x, y) coordinates
coords <- matrix(rnorm(20), ncol = 2)

#create design matrix
X <- cbind(1, coords)

#create mean for observed data to be generated
mu <- X \%*\% c(1, 2, 3)

#generate covariance matrix
V <- exp(-dist1(coords))

#generate observe data
y <- rmvnorm(mu = mu, V = V)

#find maximum likelihood estimates of covariance parameters
maxlik.cov.sp(X = X, y = y, coords = coords,
    sp.type = "exponential", range.par = 1, error.ratio = 0, 
    reml = TRUE)
}