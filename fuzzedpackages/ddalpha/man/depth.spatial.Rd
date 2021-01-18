\name{depth.spatial}
\alias{depth.spatial}
\title{
Calculate Spatial Depth
}
\description{
Calculates the spatial depth of points w.r.t. a multivariate data set.
}
\usage{
depth.spatial(x, data, mah.estimate = "moment", mah.parMcd = 0.75)
}
\arguments{
  \item{x}{
Matrix of objects (numerical vector as one object) whose depth is to be calculated; each row contains a \eqn{d}-variate point. Should have the same dimension as \code{data}.
}
  \item{data}{
Matrix of data where each row contains a \eqn{d}-variate point, w.r.t. which the depth is to be calculated.
}
  \item{mah.estimate}{ is a character string specifying which estimates to use when calculating sample covariance matrix; can be \code{"none"}, \code{"moment"} or \code{"MCD"}, determining whether traditional moment or Minimum Covariance Determinant (MCD) (see \code{\link{covMcd}}) estimates for mean and covariance are used. By default \code{"moment"} is used. With \code{"none"} the non-affine invariant version of Spatial depth is calculated
}
  \item{mah.parMcd}{
is the value of the argument \code{alpha} for the function \code{\link{covMcd}}; is used when \code{mah.estimate =} \code{"MCD"}.
}
}
\details{
Calculates spatial depth. Spatial depth (also L1-depth) is a distance-based depth exploiting the idea of spatial quantiles of Chaudhuri (1996) and Koltchinskii (1997), formulated by Vardi & Zhang (2000) and Serfling (2002).
}
\value{
Numerical vector of depths, one for each row in \code{x}; or one depth value if \code{x} is a numerical vector.
}
\references{
Chaudhuri, P. (1996). On a geometric notion of quantiles for multivariate data. \emph{Journal of the Americal Statistical Association} \bold{91} 862--872.

Koltchinskii, V.I. (1997). M-estimation, convexity and quantiles. \emph{The Annals of Statistics} \bold{25} 435--477.

Serfling, R. (2006). Depth functions in nonparametric multivariate inference. In: Liu, R., Serfling, R., Souvaine, D. (eds.), \emph{Data Depth: Robust Multivariate Analysis, Computational Geometry and Applications}, American Mathematical Society, 1--16.

Vardi, Y. and Zhang, C.H. (2000). The multivariate L1-median and associated data depth. \emph{Proceedings of the National Academy of Sciences, U.S.A.} \bold{97} 1423--1426.
}
\seealso{
\code{\link{depth.halfspace}} for calculation of the Tukey depth.

\code{\link{depth.Mahalanobis}} for calculation of Mahalanobis depth.

\code{\link{depth.projection}} for calculation of projection depth.

\code{\link{depth.simplicial}} for calculation of simplicial depth.

\code{\link{depth.simplicialVolume}} for calculation of simplicial volume depth.

\code{\link{depth.zonoid}} for calculation of zonoid depth.

\code{\link{depth.potential}} for calculation of data potential.

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
cat("Depths: ", depths, "\n")
}
\keyword{ robust }
\keyword{ multivariate }
\keyword{ nonparametric }
