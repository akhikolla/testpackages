\name{depth.space.spatial}
\alias{depth.space.spatial}
\title{
Calculate Depth Space using Spatial Depth
}
\description{
Calculates the representation of the training classes in depth space using spatial depth.
}
\usage{
depth.space.spatial(data, cardinalities, mah.estimate = "moment", mah.parMcd = 0.75)
}
\arguments{
  \item{data}{
Matrix containing training sample where each row is a \eqn{d}-dimensional object, and objects of each class are kept together so that the matrix can be thought of as containing blocks of objects representing classes.
}
  \item{cardinalities}{
Numerical vector of cardinalities of each class in \code{data}, each entry corresponds to one class.
}
  \item{mah.estimate}{ is a character string specifying which estimates to use when calculating sample covariance matrix; can be \code{"none"}, \code{"moment"} or \code{"MCD"}, determining whether traditional moment or Minimum Covariance Determinant (MCD) (see \code{\link{covMcd}}) estimates for mean and covariance are used. By default \code{"moment"} is used. With \code{"none"} the non-affine invariant version of Spatial depth is calculated
}
  \item{mah.parMcd}{
is the value of the argument \code{alpha} for the function \code{\link{covMcd}}; is used when \code{mah.estimate =} \code{"MCD"}.
}
}
\details{
The depth representation is calculated in the same way as in \code{\link{depth.spatial}}, see 'References' for more information and details.
}
\value{
Matrix of objects, each object (row) is represented via its depths (columns) w.r.t. each of the classes of the training sample; order of the classes in columns corresponds to the one in the argument \code{cardinalities}.
}
\references{
Chaudhuri, P. (1996). On a geometric notion of quantiles for multivariate data. \emph{Journal of the Americal Statistical Association} \bold{91} 862--872.

Koltchinskii, V.I. (1997). M-estimation, convexity and quantiles. \emph{The Annals of Statistics} \bold{25} 435--477.

Serfling, R. (2006). Depth functions in nonparametric multivariate inference. In: Liu, R., Serfling, R., Souvaine, D. (eds.), \emph{Data Depth: Robust Multivariate Analysis, Computational Geometry and Applications}, American Mathematical Society, 1--16.

Vardi, Y. and Zhang, C.H. (2000). The multivariate L1-median and associated data depth. \emph{Proceedings of the National Academy of Sciences, U.S.A.} \bold{97} 1423--1426.
}
\seealso{
\code{\link{ddalpha.train}} and \code{\link{ddalpha.classify}} for application, \code{\link{depth.spatial}} for calculation of spatial depth.
}
\examples{
# Generate a bivariate normal location-shift classification task
# containing 20 training objects
class1 <- mvrnorm(10, c(0,0), 
                  matrix(c(1,1,1,4), nrow = 2, ncol = 2, byrow = TRUE))
class2 <- mvrnorm(10, c(2,2), 
                  matrix(c(1,1,1,4), nrow = 2, ncol = 2, byrow = TRUE))
data <- rbind(class1, class2)
# Get depth space using spatial depth
depth.space.spatial(data, c(10, 10))

data <- getdata("hemophilia")
cardinalities = c(sum(data$gr == "normal"), sum(data$gr == "carrier"))
depth.space.spatial(data[,1:2], cardinalities)
}
\keyword{ robust }
\keyword{ multivariate }
\keyword{ nonparametric }
