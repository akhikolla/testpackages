\name{mModeAutoCovariance}
\alias{mModeAutoCovariance}
\title{The m-Mode Autocovariance Matrix}
\description{
Estimates the m-mode autocovariance matrix from an array of array-valued observations with the specified lag.
}
\usage{
mModeAutoCovariance(x, m, lag, center = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{Array of order higher than two with the last dimension corresponding to the sampling units.}
  \item{m}{The mode with respect to which the autocovariance matrix is to be computed.}
  \item{lag}{The lag with respect to which the autocovariance matrix is to be computed.}
  \item{center}{Logical, indicating whether the observations should be centered prior to computing the autocovariance matrix. Default is \code{TRUE}.}
}
\details{
The m-mode autocovariance matrix provides a higher order analogy for the ordinary autocovariance matrix of a random vector and is computed for a random tensor \eqn{X_t}{X_t} of size \eqn{p_1 \times p_2 \times \ldots \times p_r}{p_1 x p_2 x ... x p_r} as \eqn{Cov_{m \tau}(X_t) = E(X_t^{(m)} X_{t+\tau}^{(m)T})/(p_1 \ldots p_{m-1} p_{m+1} \ldots p_r)}{Cov_m,tau(X_t) = E(X_t(m) X_(t+tau)(m)^T)/(p_1 ... p_(m-1) p_(m+1) ... p_r)}, where \eqn{X_t^{(m)}}{X_t(m)} is the centered \eqn{m}{m}-flattening of \eqn{X_t}{X_t} and \eqn{\tau}{tau} is the desired \code{lag}. The algorithm computes the estimate of this based on the sample \code{x}.
}
\value{
The \code{m}-mode autocovariance matrix of \code{x} with respect to \code{lag} having the size \eqn{p_m \times p_m}{p_m x p_m}.
}
\references{
\cite{Virta, J. and Nordhausen, K., (2017), Blind source separation of tensor-valued time series, Signal Processing, 141, 204-216, \doi{10.1016/j.sigpro.2017.06.008}} 
}
\author{
Joni Virta
}

\seealso{
\code{\link{mModeCovariance}}
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

mModeAutoCovariance(S, m = 1, lag = 1)
mModeAutoCovariance(S, m = 1, lag = 4)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ utilities }
\keyword{ array }% __ONLY ONE__ keyword per line
