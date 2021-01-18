\name{ISOpure.model_optimize.cg_code.rminimize}
\alias{ISOpure.model_optimize.cg_code.rminimize}
\title{Minimize a differentiable multivariate function}
\description{This function is a conjugate-gradient search with interpolation/extrapolation by Carl 
	Edward Rasmussen.  A description of the Matlab code can be found at 
	http://learning.eng.cam.ac.uk/carl/code/minimize/ (accessed Jan. 21, 2014). This is a 
	implementation in R.}
\usage{
ISOpure.model_optimize.cg_code.rminimize(X, f, df, run_length, \ldots)
}
\arguments{
  \item{X}{The starting point is given by X which must be either a scalar or a column vector or
  	matrix, not a row matrix}
  \item{f}{The name of the function to be minimized, returning a scalar}
  \item{df}{The name of the function which returns the vector of partial derivatives of f wrt X,
	where again the partial derivatives must be in scalar or column vector/matrix form}
  \item{run_length}{Gives the length of the run: if it is positive, it gives the maximum number 
    of line searches, if negative its absolute gives the maximum allowed number of function 
    evaluations. Note, for ISOpureR, used only positive run_length.}
  \item{...}{Parameters to be passed on to the function f.}
}
\value{A list with three components:
	\item{X}{The found solution X}
	\item{fX}{A vector of function values fX indicating the progress made}
	\item{i}{The number of iterations}
	}
\details{The function returns when either its length is up, or if no further progress can be made (ie,
	we are at a (local) minimum, or so close that due to numerical problems, we cannot get any 
	closer). NOTE: If the function terminates within a few iterations, it could be an indication 
	that the function values and derivatives are not consistent (ie, there may be a bug in the 
	implementation of your "f" function).

	The Polack-Ribiere flavour of conjugate gradients is used to compute search directions, and a 
	line search using quadratic and cubic polynomial approximations and the Wolfe-Powell stopping 
	criteria is used together with the slope ratio method for guessing initial step sizes. 
	Additionally a bunch of checks are made to make sure that exploration is taking place and that
	extrapolation will not be unboundedly large.}
\author{Catalina Anghel, Francis Nguyen, Carl Edward Rasmussen}
\examples{
# Example from Carl E. Rasmussen's webpage

rosenbrock <- function(x){
	D <- length(x);
 	y <- sum(100*(x[2:D] - x[1:(D-1)]^2)^2 + (1-x[1:(D-1)])^2);
 	return(y);
 	};
drosenbrock <- function(x){
	D <- length(x);
	df <- numeric(D);
	df[1:D-1] <- -400*x[1:(D-1)]*(x[2:D]-x[1:(D-1)]^2) - 2*(1-x[1:(D-1)]);
  	df[2:D] <- df[2:D] + 200*(x[2:D]-x[1:(D-1)]^2);
  	return(df);
	};

ISOpure.model_optimize.cg_code.rminimize(c(0,0), rosenbrock, drosenbrock, 25)
#
# [[1]]
# [1] 1 1
#
# [[2]]
#  [1] 1.000000e+00 7.716094e-01 5.822402e-01 4.049274e-01 3.246633e-01
#  [6] 2.896041e-01 7.623420e-02 6.786212e-02 3.378424e-02 1.089908e-03
# [11] 1.087952e-03 8.974308e-05 1.218382e-07 6.756019e-09 3.870791e-15
# [16] 1.035408e-21 6.248025e-27 5.719242e-30 4.930381e-32
#
# [[3]]
# [1] 20
}
\keyword{optimize}
\keyword{math}