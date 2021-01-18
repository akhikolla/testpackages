\name{depth.space.halfspace}
\alias{depth.space.halfspace}
\title{
Calculate Depth Space using Halfspace Depth
}
\description{
Calculates the representation of the training classes in depth space using the halfspace depth.
}
\usage{
depth.space.halfspace(data, cardinalities, exact, method, num.directions = 1000, seed = 0)
}
\arguments{
  \item{data}{
Matrix containing training sample where each row is a \eqn{d}-dimensional object, and objects of each class are kept together so that the matrix can be thought of as containing blocks of objects representing classes.
}
  \item{cardinalities}{
Numerical vector of cardinalities of each class in \code{data}, each entry corresponds to one class.
}
  \item{exact}{
The type of the used method. The default is \code{exact=F}, which leads to approximate computation of the halfspace depth. For \code{exact=F}, \code{method="Sunif.1D"} is used by default. If \code{exact=T}, the halfspace depth is computed exactly, with \code{method="recursive"} by default.}
  \item{method}{
For \code{exact=F}, if \code{method="Sunif.1D"} (by default), the halfspace depth is computed approximately by being minimized over univariate projections (see details).

For \code{exact=T}, the halfspace depth is calculated as the minimum over all combinations of \eqn{k} points from \code{data} (see details). In this case parameter \code{method} specifies \eqn{k}, with possible values \eqn{1} for \code{method="recursive"} (by default), \eqn{d-2} for \code{method="plane"}, \eqn{d-1} for \code{method="line"}.

The name of the method may be given as well as just parameter \code{exact}, in which case the default method will be used.
}
  \item{num.directions}{
Number of random directions to be generated. As the same direction set is used for all observations, the algorithmic complexity of calculating the depth of each single point in \code{data} is logarithmic in the number of observations in \code{data}, given the number of directions, see Mozharovskyi et al. (2015), Section 2.3 for discussion.
}
  \item{seed}{
The random seed. The default value \code{seed=0} makes no changes.
}
}
\details{
The depth representation is calculated in the same way as in \code{\link{depth.halfspace}}, see References below for more information and details.
}
\value{
Matrix of objects, each object (row) is represented via its depths (columns) w.r.t. each of the classes of the training sample; order of the classes in columns corresponds to the one in the argument \code{cardinalities}.
}
\references{
Cuesta-Albertos, J.A. and Nieto-Reyes, A. (2008). The random Tukey depth. \emph{Computational Statistics and Data Analysis} \bold{52} 4979--4988.

Dyckerhoff, R. and Mozharovskyi, P. (2016). Exact computation of the halfspace depth. \emph{Computational Statistics and Data Analysis} \bold{98} 19--30.

Mozharovskyi, P., Mosler, K., and Lange, T. (2015). Classifying real-world data with the DD\eqn{\alpha}-procedure. \emph{Advances in Data Analysis and Classification} \bold{9} 287--314.

Rousseeuw, P.J. and Ruts, I. (1996). Algorithm AS 307: Bivariate location depth. \emph{Journal of the Royal Statistical Society. Series C (Applied Statistics)} \bold{45} 516--526.

Tukey, J.W. (1974). Mathematics and the picturing of data. In: \emph{Proceeding of the International Congress of Mathematicians}, Vancouver, 523--531.
}
\seealso{
\code{\link{ddalpha.train}} and \code{\link{ddalpha.classify}} for application, \code{\link{depth.halfspace}} for calculation of the Tukey depth.
}
\examples{
# Generate a bivariate normal location-shift classification task
# containing 20 training objects
class1 <- mvrnorm(10, c(0,0), 
                  matrix(c(1,1,1,4), nrow = 2, ncol = 2, byrow = TRUE))
class2 <- mvrnorm(10, c(1,1), 
                  matrix(c(1,1,1,4), nrow = 2, ncol = 2, byrow = TRUE))
data <- rbind(class1, class2)
plot(data, col = c(rep(1,10), rep(2,10)))
# Get depth space using the random Tukey depth
dhA = depth.space.halfspace(data, c(10, 10))
(dhA)

# Get depth space using default exact method - "recursive"
dhE = depth.space.halfspace(data, c(10, 10), exact = TRUE)
(dhE)

data <- getdata("hemophilia")
cardinalities = c(sum(data$gr == "normal"), sum(data$gr == "carrier"))
depth.space.halfspace(data[,1:2], cardinalities)
}
\keyword{ robust }
\keyword{ multivariate }
\keyword{ nonparametric }
