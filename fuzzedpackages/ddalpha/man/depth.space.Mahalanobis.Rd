\name{depth.space.Mahalanobis}
\alias{depth.space.Mahalanobis}
\title{
Calculate Depth Space using Mahalanobis Depth
}
\description{
Calculates the representation of the training classes in depth space using Mahalanobis depth.
}
\usage{
depth.space.Mahalanobis(data, cardinalities, mah.estimate = "moment", mah.parMcd = 0.75)
}
\arguments{
  \item{data}{
Matrix containing training sample where each row is a \eqn{d}-dimensional object, and objects of each class are kept together so that the matrix can be thought of as containing blocks of objects representing classes.
}
  \item{cardinalities}{
Numerical vector of cardinalities of each class in \code{data}, each entry corresponds to one class.
}
  \item{mah.estimate}{ is a character string specifying which estimates to use when calculating the Mahalanobis depth; can be \code{"moment"} or \code{"MCD"}, determining whether traditional moment or Minimum Covariance Determinant (MCD) (see \code{\link{covMcd}}) estimates for mean and covariance are used. By default \code{"moment"} is used.
}
  \item{mah.parMcd}{
is the value of the argument \code{alpha} for the function \code{\link{covMcd}}; is used when \code{mah.estimate =} \code{"MCD"}.
}
}
\details{
The depth representation is calculated in the same way as in \code{\link{depth.Mahalanobis}}, see 'References' for more information and details.
}
\value{
Matrix of objects, each object (row) is represented via its depths (columns) w.r.t. each of the classes of the training sample; order of the classes in columns corresponds to the one in the argument \code{cardinalities}.
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
\code{\link{ddalpha.train}} and \code{\link{ddalpha.classify}} for application, \code{\link{depth.Mahalanobis}} for calculation of Mahalanobis depth.
}
\examples{
# Generate a bivariate normal location-shift classification task
# containing 20 training objects
class1 <- mvrnorm(10, c(0,0), 
                  matrix(c(1,1,1,4), nrow = 2, ncol = 2, byrow = TRUE))
class2 <- mvrnorm(10, c(2,2), 
                  matrix(c(1,1,1,4), nrow = 2, ncol = 2, byrow = TRUE))
data <- rbind(class1, class2)
# Get depth space using Mahalanobis depth
depth.space.Mahalanobis(data, c(10, 10))
depth.space.Mahalanobis(data, c(10, 10), mah.estimate = "MCD", mah.parMcd = 0.75)

data <- getdata("hemophilia")
cardinalities = c(sum(data$gr == "normal"), sum(data$gr == "carrier"))
depth.space.Mahalanobis(data[,1:2], cardinalities)
}
\keyword{ robust }
\keyword{ multivariate }
\keyword{ nonparametric }
