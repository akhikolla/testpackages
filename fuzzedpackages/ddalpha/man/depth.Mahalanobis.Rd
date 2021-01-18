\name{depth.Mahalanobis}
\alias{depth.Mahalanobis}
\title{
Calculate Mahalanobis Depth
}
\description{
Calculates the Mahalanobis depth of points w.r.t. a multivariate data set.
}
\usage{
depth.Mahalanobis(x, data, mah.estimate = "moment", mah.parMcd = 0.75)
}
\arguments{
  \item{x}{
Matrix of objects (numerical vector as one object) whose depth is to be calculated; each row contains a \eqn{d}-variate point. Should have the same dimension as \code{data}.
}
  \item{data}{
Matrix of data where each row contains a \eqn{d}-variate point, w.r.t. which the depth is to be calculated.
}
  \item{mah.estimate}{ is a character string specifying which estimates to use when calculating the Mahalanobis depth; can be \code{"moment"} or \code{"MCD"}, determining whether traditional moment or Minimum Covariance Determinant (MCD) (see \code{\link{covMcd}}) estimates for mean and covariance are used. By default \code{"moment"} is used.
}
  \item{mah.parMcd}{
is the value of the argument \code{alpha} for the function \code{\link{covMcd}}; is used when \code{mah.estimate =} \code{"MCD"}.
}
}
\details{
Calculates Mahalanobis depth. Mahalanobis depth is based on an outlyingness measure (Zuo & Serfling, 2000), \emph{viz.} the Mahalanobis distance between the given point and the center of the data (Mahalanobis, 1936).

 \emph{Moment estimates} may be used i.e. traditional \emph{mean} and \emph{covariance matrix}, the corresponding depth may be sensitive to
outliers. A more robust depth is obtained with \emph{minimum volume ellipsoid} (MVE) or \emph{minimum
covariance determinant} (MCD) estimators, see Rousseeuw & Leroy (1987) and Lopuhaa &
Rousseeuw (1991).
}
\value{
Numerical vector of depths, one for each row in \code{x}; or one depth value if \code{x} is a numerical vector.
}
\references{
Mahalanobis, P. (1936). On the generalized distance in statistics. \emph{Proceedings of the National
Academy India} \bold{12} 49--55.

Liu, R.Y. (1992). Data depth and multivariate rank tests. In: Dodge, Y. (ed.), \emph{L1-Statistics and Related Methods}, North-Holland (Amsterdam), 279--294.

Lopuhaa, H.P. and Rousseeuw, P.J. (1991). Breakdown points of affine equivariant estimators of multivariate location and covariance matrices. \emph{The Annals of Statistics} \bold{19} 229--248.

Rousseeuw, P.J. and Leroy, A.M. (1987). Robust Regression and Outlier Detection. John Wiley & Sons (New York).

Zuo, Y.J. and Serfling, R. (2000). General notions of statistical depth function. \emph{The Annals of Statistics} \bold{28} 461--482.
}
\seealso{
\code{\link{depth.halfspace}} for calculation of the Tukey depth.

\code{\link{depth.projection}} for calculation of projection depth.

\code{\link{depth.simplicial}} for calculation of simplicial depth.

\code{\link{depth.simplicialVolume}} for calculation of simplicial volume depth.

\code{\link{depth.spatial}} for calculation of spatial depth.

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
                
depths <- depth.Mahalanobis(x, data)
cat("Depths moment: ", depths, "\n")
depths <- depth.Mahalanobis(x, data, mah.estimate = "MCD", mah.parMcd = 0.75)
cat("Depths MCD: ", depths, "\n")
}
\keyword{ robust }
\keyword{ multivariate }
\keyword{ nonparametric }
