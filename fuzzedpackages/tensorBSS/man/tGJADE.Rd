\name{tgJADE}
\alias{tgJADE}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{gJADE for Tensor-Valued Time Series}
\description{
Computes the tensorial gJADE for time series where at each time point a tensor of order \eqn{r} is observed. 
}
\usage{
tgJADE(x, lags = 0:12, maxiter = 100, eps = 1e-06)
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
The assumption is that the elements of \eqn{S} are mutually independent, centered and weakly stationary time series and are mixed from each mode
\eqn{m} by the mixing matrix \eqn{A_m}, \eqn{m = 1, \ldots, r}{m= 1, ..., r}, yielding the observed time series \eqn{X}. In R the sample of \eqn{X} is saved as an \code{\link{array}} of dimensions
\eqn{p_1, p_2, \ldots, p_r, T}{p_1, p_2, ..., p_r, T}.

\code{tgJADE} recovers then based on \code{x} the underlying independent time series \eqn{S} by estimating the \eqn{r} unmixing matrices 
 \eqn{W_1, \ldots, W_r}{W_1, ..., W_r} using the lagged fourth joint moments specified by \code{lags}. This reliance on higher order moments makes the method especially suited for stochastic volatility models.

If \code{x} is a matrix, that is, \eqn{r = 1}{r = 1}, the method reduces to gJADE and the function calls \code{\link[tsBSS]{gJADE}}.

If \code{lags = 0} the method reduces to \code{\link{tJADE}}.
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
\code{\link[tsBSS]{gJADE}}, \code{\link[JADE]{rjd}}, \code{\link{tJADE}}
}
\examples{
library("stochvol")
n <- 1000
S <- t(cbind(svsim(n, mu = -10, phi = 0.98, sigma = 0.2, nu = Inf)$y,
             svsim(n, mu = -5, phi = -0.98, sigma = 0.2, nu = 10)$y,
             svsim(n, mu = -10, phi = 0.70, sigma = 0.7, nu = Inf)$y,
             svsim(n, mu = -5, phi = -0.70, sigma = 0.7, nu = 10)$y,
             svsim(n, mu = -9, phi = 0.20, sigma = 0.01, nu = Inf)$y,
             svsim(n, mu = -9, phi = -0.20, sigma = 0.01, nu = 10)$y))
dim(S) <- c(3, 2, n)

A1 <- matrix(rnorm(9), 3, 3)
A2 <- matrix(rnorm(4), 2, 2)

X <- tensorTransform(S, A1, 1)
X <- tensorTransform(X, A2, 2)

tgjade <- tgJADE(X)

MD(tgjade$W[[1]], A1)
MD(tgjade$W[[2]], A2) 
tMD(tgjade$W, list(A1, A2))

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{multivariate }
\keyword{array }
\keyword{ts }% __ONLY ONE__ keyword per line


