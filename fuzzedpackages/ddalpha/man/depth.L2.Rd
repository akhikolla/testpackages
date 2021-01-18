\name{depth.L2}
\alias{depth.L2}
\title{
Calculate L2-Depth
}
\description{
Calculates the L2-depth of points w.r.t. a multivariate data set.
}
\usage{
depth.L2(x, data, mah.estimate = "moment", mah.parMcd = 0.75)
}
\arguments{
  \item{x}{
Matrix of objects (numerical vector as one object) whose depth is to be calculated; each row contains a \eqn{d}-variate point. Should have the same dimension as \code{data}.
}
  \item{data}{
Matrix of data where each row contains a \eqn{d}-variate point, w.r.t. which the depth is to be calculated.
}
  \item{mah.estimate}{ is a character string specifying which estimates to use when calculating sample covariance matrix; can be \code{"none"}, \code{"moment"} or \code{"MCD"}, determining whether traditional moment or Minimum Covariance Determinant (MCD) (see \code{\link{covMcd}}) estimates for mean and covariance are used. By default \code{"moment"} is used. With \code{"none"} the non-affine invariant version of the L2-depth is calculated
}
  \item{mah.parMcd}{
is the value of the argument \code{alpha} for the function \code{\link{covMcd}}; is used when \code{mah.estimate =} \code{"MCD"}.
}
}
\details{
Calculates L2-depth (Mosler, 2013). L2-depth is based on the oultyingness distance calculated as the average L2-distance from (a row of) \code{x} to each point in \code{data}.
}
\value{
Numerical vector of depths, one for each row in \code{x}; or one depth value if \code{x} is a numerical vector.
}
\references{
Mosler, K. (2013). Depth statistics. In: Becker, C., Fried, R. and Kuhnt, S. (eds), \emph{Robustness and Complex Data Structures: Festschrift in Honour of Ursula Gather}, Springer-Verlag (Berlin, Heidelberg), 17--34.
}
\seealso{
\code{\link{depth.halfspace}} for calculation of the Tukey depth.

\code{\link{depth.Mahalanobis}} for calculation of Mahalanobis depth.

\code{\link{depth.projection}} for calculation of projection depth.

\code{\link{depth.qhpeeling}} for calculation of convex hull peeling depth.

\code{\link{depth.simplicial}} for calculation of simplicial depth.

\code{\link{depth.simplicialVolume}} for calculation of simplicial volume depth.

\code{\link{depth.spatial}} for calculation of spatial depth.

\code{\link{depth.potential}} for calculation of data potential.

\code{\link{depth.zonoid}} for calculation of zonoid depth.

}
\examples{
# 5-dimensional normal distribution
data <- mvrnorm(1000, rep(0, 5), 
                matrix(c(1, 0, 0, 0, 0, 
                         0, 2, 0, 0, 0, 
                         0, 0, 3, 0, 0, 
                         0, 0, 0, 2, 0, 
                         0, 0, 0, 0, 1),
                nrow = 5))
x <- mvrnorm(10, rep(1, 5), 
             matrix(c(1, 0, 0, 0, 0, 
                      0, 1, 0, 0, 0, 
                      0, 0, 1, 0, 0, 
                      0, 0, 0, 1, 0, 
                      0, 0, 0, 0, 1),
             nrow = 5))
                
depths <- depth.spatial(x, data)
cat("Depths:", depths, "\n")
}
\keyword{ robust }
\keyword{ multivariate }
\keyword{ nonparametric }
