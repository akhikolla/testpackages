\name{tensorStandardize}
\alias{tensorStandardize}

\title{Standardize an Observation Array}
\description{
Standardizes an array of array-valued observations simultaneously from each mode. The method can be seen as a higher-order analogy for the regular multivariate standardization of random vectors.
}
\usage{
tensorStandardize(x, location = NULL, scatter = NULL)
}
\arguments{
  \item{x}{Array of an order higher than two with the last dimension corresponding to the sampling units.}
  \item{location}{The location to be used in the standardizing. Either \code{NULL}, defaulting to the mean array, or a user-specified \eqn{p_1 \times p_2 \times \ldots \times p_r}{p_1 x p_2 x ... x p_r}-dimensional array.}
  \item{scatter}{The scatter matrices to be used in the standardizing. Either \code{NULL}, defaulting to the m-mode covariance matrices, or a user-specified list of length \code{r} of \eqn{p_1 \times p_1, \ldots , p_r \times p_r}{p_1 x p_1 ... p_r x p_r}-dimensional symmetric positive definite matrices.}
}
\details{
The algorithm first centers the \eqn{n}{n} observed tensors \eqn{X_i}{Xi} using \code{location} (either the sample mean, or a user-specified location). Then, if \code{scatter = NULL}, it estimates the \eqn{m}{m}th mode covariance matrix \eqn{Cov_m(X) = E(X^{(m)} X^{(m)T})/(p_1 \ldots p_{m-1} p_{m+1} \ldots p_r)}{Cov_m(X) = E(X(m) X(m)^T)/(p1 ... p(m-1) p(m+1) ... pr)}, where \eqn{X^{(m)}}{X(m)} is the centered \eqn{m}{m}-flattening of \eqn{X}{X}, for each mode, and transforms the observations with the inverse square roots of the covariance matrices from the corresponding modes. If, instead, the user has specified a non-\code{NULL} value for \code{scatter}, the inverse square roots of those matrices are used to transform the centered data.
}
\value{
A list containing the following components: 
\item{x }{Array of the same size as \code{x} containing the standardized observations. The used location and scatters are returned as attributes \code{"location"} and \code{"scatter"}.}
\item{S }{List containing inverse square roots of the covariance matrices of different modes.}
}
\author{
Joni Virta
}
\examples{
# Generate sample data.
n <- 100
x <- t(cbind(rnorm(n, mean = 0),
             rnorm(n, mean = 1),
             rnorm(n, mean = 2),
             rnorm(n, mean = 3),
             rnorm(n, mean = 4),
             rnorm(n, mean = 5)))
             
dim(x) <- c(3, 2, n)

# Standardize
z <- tensorStandardize(x)$x

# The m-mode covariance matrices of the standardized tensors
mModeCovariance(z, 1)
mModeCovariance(z, 2)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ utilities }
\keyword{ array }
