\name{depth.space.zonoid}
\alias{depth.space.zonoid}
\title{
Calculate Depth Space using Zonoid Depth
}
\description{
Calculates the representation of the training classes in depth space using zonoid depth.
}
\usage{
depth.space.zonoid(data, cardinalities, seed = 0)
}
\arguments{
  \item{data}{
Matrix containing training sample where each row is a \eqn{d}-dimensional object, and objects of each class are kept together so that the matrix can be thought of as containing blocks of objects representing classes.
}
  \item{cardinalities}{
Numerical vector of cardinalities of each class in \code{data}, each entry corresponds to one class.
}
  \item{seed}{
the random seed. The default value \code{seed=0} makes no changes.
}
}
\details{
The depth representation is calculated in the same way as in \code{\link{depth.zonoid}}, see 'References' for more information and details.
}
\value{
Matrix of objects, each object (row) is represented via its depths (columns) w.r.t. each of the classes of the training sample; order of the classes in columns corresponds to the one in the argument \code{cardinalities}.
}
\references{
Dyckerhoff, R., Koshevoy, G., and Mosler, K. (1996). Zonoid data depth: theory and computation. In: Prat A. (ed), \emph{COMPSTAT 1996. Proceedings in computational statistics}, Physica-Verlag (Heidelberg), 235--240.

Koshevoy, G. and Mosler, K. (1997). Zonoid trimming for multivariate distributions \emph{Annals of Statistics} \bold{25} 1998--2017.

Mosler, K. (2002). \emph{Multivariate dispersion, central regions and depth: the lift zonoid approach} Springer (New York).
}
\seealso{
\code{\link{ddalpha.train}} and \code{\link{ddalpha.classify}} for application, \code{\link{depth.zonoid}} for calculation of zonoid depth.
}
\examples{
# Generate a bivariate normal location-shift classification task
# containing 20 training objects
class1 <- mvrnorm(10, c(0,0), 
                  matrix(c(1,1,1,4), nrow = 2, ncol = 2, byrow = TRUE))
class2 <- mvrnorm(10, c(2,2), 
                  matrix(c(1,1,1,4), nrow = 2, ncol = 2, byrow = TRUE))
data <- rbind(class1, class2)
# Get depth space using zonoid depth
depth.space.zonoid(data, c(10, 10))

data <- getdata("hemophilia")
cardinalities = c(sum(data$gr == "normal"), sum(data$gr == "carrier"))
depth.space.zonoid(data[,1:2], cardinalities)

}
\keyword{ robust }
\keyword{ multivariate }
\keyword{ nonparametric }
