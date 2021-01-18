\name{depth.simplicialVolume}
\alias{depth.simplicialVolume}
\title{
Calculate Simplicial Volume Depth
}
\description{
Calculates the simpicial volume depth of points w.r.t. a multivariate data set.
}
\usage{
depth.simplicialVolume(x, data, exact = F, k = 0.05, mah.estimate = "moment", 
                       mah.parMcd = 0.75, seed = 0)
}
\arguments{
  \item{x}{
Matrix of objects (numerical vector as one object) whose depth is to be calculated; each row contains a \eqn{d}-variate point. Should have the same dimension as \code{data}.
}
  \item{data}{
Matrix of data where each row contains a \eqn{d}-variate point, w.r.t. which the depth is to be calculated.
}
  \item{exact}{
\code{exact=F} (by default) implies the approximative algorithm, considering \code{k} simplices, \code{exact=T} implies the exact algorithm.
}
  \item{k}{
Number (\eqn{k>1}) or portion (if \eqn{0<k<1}) of simplices that are considered if \code{exact=F}. If \eqn{k>1}, then the algorithmic complexity is polynomial in \eqn{d} but is independent of the number of observations in \code{data}, given \eqn{k}. If \eqn{0<k<1}, then the algorithmic complexity is exponential in the number of observations in \code{data}, but the calculation precision stays approximately the same.
}
  \item{mah.estimate}{
A character string specifying affine-invariance adjustment; can be \code{"none"}, \code{"moment"} or \code{"MCD"}, determining whether no affine-invariance adjustemt or moment or Minimum Covariance Determinant (MCD) (see \code{\link{covMcd}}) estimates of the covariance are used. By default \code{"moment"} is used.
}
  \item{mah.parMcd}{
The value of the argument \code{alpha} for the function \code{\link{covMcd}}; is used when \code{mah.estimate =} \code{"MCD"}.
}
  \item{seed}{
The random seed. The default value \code{seed=0} makes no changes.
}
}
\details{
Calculates Oja depth (also: Simplicial volume depth). 
At first the Oja outlyingness function \code{O(x,data)} is calculated as the average of the volumes of simplices built on \eqn{d} data points and the measurement point \code{x} (Oja, 1983).

Zuo and Serfling (2000) proposed Oja depth based on the Oja outlyingness function as \code{1/(1 + O(x,data)/S)}, where S is a square root of the determinant of \code{cov(data)}, which makes the depth function affine-invariant.
}
\value{
Numerical vector of depths, one for each row in \code{x}; or one depth value if \code{x} is a numerical vector.
}
\references{
Oja, H. (1983). Descriptive statistics for multivariate distributions. \emph{Statistics & Probability Letters} \bold{1} 327--332.

Zuo, Y.J. and Serfling, R. (2000). General notions of statistical depth function. \emph{The Annals of Statistics} \bold{28} 461--482.
}
\seealso{
\code{\link{depth.halfspace}} for calculation of the Tukey depth.

\code{\link{depth.Mahalanobis}} for calculation of Mahalanobis depth.

\code{\link{depth.projection}} for calculation of projection depth.

\code{\link{depth.simplicial}} for calculation of simplicial depth.

\code{\link{depth.spatial}} for calculation of spatial depth.

\code{\link{depth.zonoid}} for calculation of zonoid depth.

\code{\link{depth.potential}} for calculation of data potential.

}
\examples{
# 3-dimensional normal distribution
data <- mvrnorm(20, rep(0, 3), 
                matrix(c(1, 0, 0,
                         0, 2, 0,
                         0, 0, 1),
                       nrow = 3))
x <- mvrnorm(10, rep(1, 3), 
             matrix(c(1, 0, 0,
                      0, 1, 0,
                      0, 0, 1),
                    nrow = 3))

#exact
depths <- depth.simplicialVolume(x, data, exact = TRUE)
cat("Depths: ", depths, "\n")

#approximative
depths <- depth.simplicialVolume(x, data, exact = FALSE, k = 0.2)
cat("Depths: ", depths, "\n")
}
\keyword{ robust }
\keyword{ multivariate }
\keyword{ nonparametric }
