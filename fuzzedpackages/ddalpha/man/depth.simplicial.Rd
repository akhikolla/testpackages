\name{depth.simplicial}
\alias{depth.simplicial}
\title{
Calculate Simplicial Depth
}
\description{
Calculates the simplicial depth of points w.r.t. a multivariate data set.
}
\usage{
depth.simplicial(x, data, exact = F, k = 0.05, seed = 0)
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
  \item{seed}{
the random seed. The default value \code{seed=0} makes no changes.
}
}
\details{
Calculates simplicial depth. Simplicial depth is counted as a probability that a point lies in a simplex, built on \eqn{d+1} data points.
}
\value{
Numerical vector of depths, one for each row in \code{x}; or one depth value if \code{x} is a numerical vector.
}
\references{
Chaudhuri, P. (1996). On a geometric notion of quantiles for multivariate data. \emph{Journal of the American Statistical Association} \bold{91} 862--872.

Liu, R. Y. (1990). On a notion of data depth based on random simplices. \emph{The Annals of Statistics} \bold{18} 405--414.

Rousseeuw, P.J. and Ruts, I. (1996). Algorithm AS 307: Bivariate location depth. \emph{Journal of the Royal Statistical Society. Seriec C (Applied Statistics)} \bold{45} 516--526.
}
\seealso{
\code{\link{depth.halfspace}} for calculation of the Tukey depth.

\code{\link{depth.Mahalanobis}} for calculation of Mahalanobis depth.

\code{\link{depth.projection}} for calculation of projection depth.

\code{\link{depth.simplicialVolume}} for calculation of simplicial volume depth.

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
depths <- depth.simplicial(x, data, exact = TRUE)
cat("Depths: ", depths, "\n")

#approximative
depths <- depth.simplicial(x, data, exact = FALSE, k = 0.2)
cat("Depths: ", depths, "\n")
}
\keyword{ robust }
\keyword{ multivariate }
\keyword{ nonparametric }
