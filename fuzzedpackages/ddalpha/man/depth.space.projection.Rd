\name{depth.space.projection}
\alias{depth.space.projection}
\title{
Calculate Depth Space using Projection Depth
}
\description{
Calculates the representation of the training classes in depth space using projection depth.
}
\usage{
depth.space.projection(data, cardinalities, 
                       method = "random", num.directions = 1000, seed = 0)
}
\arguments{
  \item{data}{
Matrix containing training sample where each row is a \eqn{d}-dimensional object, and objects of each class are kept together so that the matrix can be thought of as containing blocks of objects representing classes.
}
  \item{cardinalities}{
Numerical vector of cardinalities of each class in \code{data}, each entry corresponds to one class.
}
  \item{method}{
  to be used in calculations.
  
  \code{"random"} Here the depth is determined as the minimum univariate depth of the data projected on lines in several directions. The directions are distributed uniformly on the \eqn{(d-1)}-sphere; the same direction set is used for all points.
  
  \code{"linearize"} The Nelder-Mead method for function minimization, taken from Olsson, Journal of Quality Technology, 1974, 6, 56. R-codes of this function were written by Subhajit Dutta.
  }
  \item{num.directions}{
Number of random directions to be generated for \code{method = "random"}. With the growth of n the complexity grows linearly for the same number of directions. 
}
  \item{seed}{
the random seed. The default value \code{seed=0} makes no changes.
}
}
\details{
The depth representation is calculated in the same way as in \code{\link{depth.projection}}, see 'References' for more information and details.
}
\value{
Matrix of objects, each object (row) is represented via its depths (columns) w.r.t. each of the classes of the training sample; order of the classes in columns corresponds to the one in the argument \code{cardinalities}.
}
\references{
Donoho, D.L. (1982). \emph{Breakdown properties of multivariate location estimators}. Ph.D. qualifying paper. Department of Statistics, Harvard University.

Liu, R.Y. (1992). Data depth and multivariate rank tests. In: Dodge, Y. (ed.), L1-Statistics and Related Methods, North-Holland (Amsterdam), 279--294.

Liu, X. and Zuo, Y. (2014). Computing projection depth and its associated estimators. \emph{Statistics and Computing} \bold{24} 51--63.

Stahel, W.A. (1981). \emph{Robust estimation: infinitesimal optimality and covariance matrix estimators}. Ph.D. thesis (in German). Eidgenossische Technische Hochschule Zurich.

Zuo, Y.J. and Lai, S.Y. (2011). Exact computation of bivariate projection depth and the Stahel-Donoho estimator. \emph{Computational Statistics and Data Analysis} \bold{55} 1173--1179.
}
\seealso{
\code{\link{ddalpha.train}} and \code{\link{ddalpha.classify}} for application, \code{\link{depth.projection}} for calculation of projection depth.
}
\examples{
# Generate a bivariate normal location-shift classification task
# containing 20 training objects
class1 <- mvrnorm(10, c(0,0), 
                  matrix(c(1,1,1,4), nrow = 2, ncol = 2, byrow = TRUE))
class2 <- mvrnorm(10, c(2,2), 
                  matrix(c(1,1,1,4), nrow = 2, ncol = 2, byrow = TRUE))
data <- rbind(class1, class2)
# Get depth space using projection depth
depth.space.projection(data, c(10, 10), method = "random", num.directions = 1000)
depth.space.projection(data, c(10, 10), method = "linearize")

data <- getdata("hemophilia")
cardinalities = c(sum(data$gr == "normal"), sum(data$gr == "carrier"))
depth.space.projection(data[,1:2], cardinalities)
}
\keyword{ robust }
\keyword{ multivariate }
\keyword{ nonparametric }
