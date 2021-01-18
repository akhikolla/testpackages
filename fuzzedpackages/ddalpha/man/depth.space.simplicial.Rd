\name{depth.space.simplicial}
\alias{depth.space.simplicial}
\title{
Calculate Depth Space using Simplicial Depth
}
\description{
Calculates the representation of the training classes in depth space using simplicial depth.
}
\usage{
depth.space.simplicial(data, cardinalities, exact = F, k = 0.05, seed = 0)
}
\arguments{
  \item{data}{
Matrix containing training sample where each row is a \eqn{d}-dimensional object, and objects of each class are kept together so that the matrix can be thought of as containing blocks of objects representing classes.
}
  \item{cardinalities}{
Numerical vector of cardinalities of each class in \code{data}, each entry corresponds to one class.
}
  \item{exact}{
\code{exact=F} (by default) implies the approximative algorithm, considering \code{k} simplices, \code{exact=T} implies the exact algorithm.
}
  \item{k}{
Number (\eqn{k>1}) or portion (if \eqn{0<k<1}) of simplices that are considered if \code{exact=F}. If \eqn{k>1}, then the algorithmic complexity is polynomial in \eqn{d} but is independent of the number of observations in \code{data}, given \eqn{k}. If \eqn{0<k<1}, then the algorithmic complexity is exponential in the number of observations in \code{data}, but the calculation precision stays approximately the same.
}
  \item{seed}{
The random seed. The default value \code{seed=0} makes no changes.
}
}
\details{
The depth representation is calculated in the same way as in \code{\link{depth.simplicial}}, see 'References' for more information and details.
}
\value{
Matrix of objects, each object (row) is represented via its depths (columns) w.r.t. each of the classes of the training sample; order of the classes in columns corresponds to the one in the argument \code{cardinalities}.
}
\references{
Chaudhuri, P. (1996). On a geometric notion of quantiles for multivariate data. \emph{Journal of the American Statistical Association} \bold{91} 862--872.

Liu, R. Y. (1990). On a notion of data depth based on random simplices. \emph{The Annals of Statistics} \bold{18} 405--414.

Rousseeuw, P.J. and Ruts, I. (1996). Algorithm AS 307: Bivariate location depth. \emph{Journal of the Royal Statistical Society. Seriec C (Applied Statistics)} \bold{45} 516--526.
}
\seealso{
\code{\link{ddalpha.train}} and \code{\link{ddalpha.classify}} for application, \code{\link{depth.simplicial}} for calculation of simplicial depth.
}
\examples{
# Generate a bivariate normal location-shift classification task
# containing 20 training objects
class1 <- mvrnorm(10, c(0,0), 
                  matrix(c(1,1,1,4), nrow = 2, ncol = 2, byrow = TRUE))
class2 <- mvrnorm(10, c(1,1), 
                  matrix(c(1,1,1,4), nrow = 2, ncol = 2, byrow = TRUE))
data <- rbind(class1, class2)
# Get depth space using simplicial depth
depth.space.simplicial(data, c(10, 10))

data <- getdata("hemophilia")
cardinalities = c(sum(data$gr == "normal"), sum(data$gr == "carrier"))
depth.space.simplicial(data[,1:2], cardinalities)
}
\keyword{ robust }
\keyword{ multivariate }
\keyword{ nonparametric }
