\name{tSOBI}
\alias{tSOBI}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{SOBI for Tensor-Valued Time Series
}
\description{
Computes the tensorial SOBI for time series where at each time point a tensor of order \eqn{r} is observed. 
}
\usage{
tSOBI(x, lags = 1:12, maxiter = 100, eps = 1e-06)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{Numeric array of an order at least two. It is assumed that the last dimension corresponds to the time.}
  \item{lags}{Vector of integers. Defines the lags used for the computations of the autocovariances.}
  \item{maxiter}{Maximum number of iterations. Passed on to \code{\link[JADE]{rjd}}.}
  \item{eps}{Convergence tolerance. Passed on to \code{\link[JADE]{rjd}}.}
}

\details{
It is assumed that \eqn{S} is a tensor (array) of size \eqn{p_1 \times p_2 \times \ldots \times p_r}{p_1 x p_2 x ... x p_r} measured at time points \eqn{1, \ldots, T}{1, ..., T}.
The assumption is that the elements of \eqn{S} are uncorrelated, centered and weakly stationary time series and are mixed from each mode
\eqn{m} by the mixing matrix \eqn{A_m}, \eqn{m = 1, \ldots, r}{m= 1, ..., r}, yielding the observed time series \eqn{X}. In R the sample of \eqn{X} is saved as an \code{\link{array}} of dimensions
\eqn{p_1, p_2, \ldots, p_r, T}{p_1, p_2, ..., p_r, T}.

\code{tSOBI} recovers then based on \code{x} the underlying uncorrelated time series \eqn{S} by estimating the \eqn{r} unmixing matrices 
 \eqn{W_1, \ldots, W_r}{W_1, ..., W_r} using the lagged joint autocovariances specified by \code{lags}.

If \code{x} is a matrix, that is, \eqn{r = 1}{r = 1}, the method reduces to SOBI and the function calls \code{\link[JADE]{SOBI}}.
}

\value{
A list with class 'tbss', inheriting from class 'bss', containing the following components: 
\item{S }{Array of the same size as x containing the estimated uncorrelated sources.}
\item{W }{List containing all the unmixing matrices}
\item{Xmu }{The data location.}
\item{datatype}{Character string with value "ts". Relevant for \code{\link{plot.tbss}}.}
}
\references{
\cite{Virta, J. and Nordhausen, K., (2017), Blind source separation of tensor-valued time series. Signal Processing 141, 204-216, \doi{10.1016/j.sigpro.2017.06.008}} 
}
\author{
Joni Virta
}



\seealso{
 \code{\link[JADE]{SOBI}}, \code{\link[JADE]{rjd}}
}
\examples{
n <- 1000
S <- t(cbind(as.vector(arima.sim(n = n, list(ar = 0.9))),
             as.vector(arima.sim(n = n, list(ar = -0.9))),
             as.vector(arima.sim(n = n, list(ma = c(0.5, -0.5)))),
             as.vector(arima.sim(n = n, list(ar = c(-0.5, -0.3)))),
             as.vector(arima.sim(n = n, list(ar = c(0.5, -0.3, 0.1, -0.1), ma=c(0.7, -0.3)))),
             as.vector(arima.sim(n = n, list(ar = c(-0.7, 0.1), ma = c(0.9, 0.3, 0.1, -0.1))))))
dim(S) <- c(3, 2, n)

A1 <- matrix(rnorm(9), 3, 3)
A2 <- matrix(rnorm(4), 2, 2)

X <- tensorTransform(S, A1, 1)
X <- tensorTransform(X, A2, 2)

tsobi <- tSOBI(X)

MD(tsobi$W[[1]], A1)
MD(tsobi$W[[2]], A2) 
tMD(tsobi$W, list(A1, A2))

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{multivariate }
\keyword{array }
\keyword{ts }% __ONLY ONE__ keyword per line
