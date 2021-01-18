\name{mModeCovariance}
\alias{mModeCovariance}
\title{The m-Mode Covariance Matrix}
\description{
Estimates the m-mode covariance matrix from an array of array-valued observations.
}
\usage{
mModeCovariance(x, m, center = TRUE)
}
\arguments{
  \item{x}{Array of order higher than two with the last dimension corresponding to the sampling units.}
  \item{m}{The mode with respect to which the covariance matrix is to be computed.}
  \item{center}{Logical, indicating whether the observations should be centered prior to computing the covariance matrix. Default is \code{TRUE}.}
}
\details{
The m-mode covariance matrix provides a higher order analogy for the ordinary covariance matrix of a random vector and is computed for a random tensor \eqn{X}{X} of size \eqn{p_1 \times p_2 \times \ldots \times p_r}{p_1 x p_2 x ... x p_r} as \eqn{Cov_m(X) = E(X^{(m)} X^{(m)T})/(p_1 \ldots p_{m-1} p_{m+1} \ldots p_r)}{Cov_m(X) = E(X(m) X(m)^T)/(p_1 ... p_(m-1) p_(m+1) ... p_r)}, where \eqn{X^{(m)}}{X(m)} is the centered \eqn{m}{m}-flattening of \eqn{X}{X}. The algorithm computes the estimate of this based on the sample \code{x}.
}
\value{
The \code{m}-mode covariance matrix of \code{x} having the size \eqn{p_m \times p_m}{p_m x p_m}.
}
\references{
\cite{Virta, J., Li, B., Nordhausen, K. and Oja, H., (2017), Independent component analysis for tensor-valued data, Journal of Multivariate Analysis, \doi{10.1016/j.jmva.2017.09.008}}
}
\author{
Joni Virta
}
\seealso{
\code{\link{mModeAutoCovariance}}
}
\examples{
## Generate sample data.
n <- 100
x <- t(cbind(rnorm(n, mean = 0),
             rnorm(n, mean = 1),
             rnorm(n, mean = 2),
             rnorm(n, mean = 3),
             rnorm(n, mean = 4),
             rnorm(n, mean = 5)))
             
dim(x) <- c(3, 2, n)

# The m-mode covariance matrices of the first and second modes
mModeCovariance(x, 1)
mModeCovariance(x, 2)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ utilities }
\keyword{ array }% __ONLY ONE__ keyword per line
