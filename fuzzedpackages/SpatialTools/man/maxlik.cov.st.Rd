\name{maxlik.cov.st} 

\alias{maxlik.cov.st}

\title{
	Determines maximum likelihood estimates of covariance parameters
} 

\description{
	Estimates covariance parameters of spatio-temporal covariance functions using maximum likelihood or restricted maximum likelihood.  See \code{cov.st} for more details of covariance functions to be estimated.  The covariance function is reparameterized slightly to speedup computation.  Specifically, the variance parameter for the hidden process, sp.par[1], is profiled out and the error.var parameter is parameterized as sp.par[1] * error.ratio.
} 

\usage{
	maxlik.cov.st(X, y, coords, time, sp.type = "exponential", 
	    range.par = stop("specify range.par argument"), 
	    error.ratio = stop("specify error.ratio argument"), 
	    smoothness = 0.5, t.type = "ar1", t.par = .5, D = NULL, T = NULL, 
	    reml = TRUE, lower = NULL, upper = NULL, control = list(trace = TRUE), 
	    optimizer="nlminb")
} 

\arguments{
	\item{X}{A numeric matrix of size \eqn{n \times k} containing the design matrix of the data locations.} 
	\item{y}{A vector of length \eqn{n} containing the observed responses.} 
	\item{coords}{A numeric matrix of size \eqn{n \times d} containing the locations of the observed responses.} 
	\item{time}{A numeric vector of length n containing the time at which the responses were observed.} 
	\item{sp.type}{A character vector specifying the spatial covariance type.  Valid types are currently exponential, gaussian, matern, and spherical.} 
	\item{range.par}{An initial guess for the spatial dependence parameter.}
	\item{error.ratio}{A value non-negative value indicating the ratio \code{error.var/sp.par[1]}.} 
	\item{smoothness}{A positive number indicating the variance of the error term.} 
	\item{t.type}{A character vector indicating the spatial covariance type.  Only \code{ar1} is currently available.} 
	\item{t.par}{A value specifying the temporal dependence parameter of the ar1 process.} 
	\item{D}{The Euclidean distance matrix for the coords matrix.  Must be of size \eqn{n \times n}.}
	\item{T}{The Euclidean distance matrix for the time matrix.  Must be of size \eqn{n \times n}.}
	\item{reml}{A boolean value indicating whether restricted maximum likelihood estimation should be used.  Defaults to TRUE.}
	\item{lower}{A vector giving lower bounds for the covariance parameters \code{sp.par[2]}, \code{error.ratio}, and \code{smoothness} (when the model is matern).  Order matters!  If not given defaults to a lower bound of .001 for \code{sp.par[2]}, 0 for \code{error.ratio}, and .001 for \code{smoothness}.}
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
	When doing the numerical optimization, the covariance function is reparameterized slightly to speedup computation.  
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
#Generate locations and observed times
coords <- matrix(rnorm(40), ncol = 2)
time <- rep(1:2, each = 10)

#Calculate distance matrix for time vector
T <- dist1(matrix(time))

#create design matrix
X <- cbind(1, coords)

#create mean for observed data to be generated
mu <- X \%*\% c(1, 2, 3)

#generate covariance matrix for spatio-temporal data
V <- exp(-dist1(coords)) * .25^T

#generate observe data
y <- rmvnorm(mu = mu, V = V)

maxlik.cov.st(X = X, y = y, coords = coords, time = time,
    sp.type = "exponential", range.par = 1, error.ratio = 0, 
    t.type = "ar1", t.par = .5, reml = TRUE)

}