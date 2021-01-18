\name{depth.space.potential}
\alias{depth.space.potential}
\title{
Calculate Potential Space
}
\description{
Calculates the representation of the training classes in potential space.
}
\usage{
depth.space.potential(data, cardinalities, pretransform = "NMom", 
            kernel = "GKernel", kernel.bandwidth = NULL, mah.parMcd = 0.75)
}
\arguments{
  \item{data}{
Matrix containing training sample where each row is a \eqn{d}-dimensional object, and objects of each class are kept together so that the matrix can be thought of as containing blocks of objects representing classes.
}
  \item{cardinalities}{
Numerical vector of cardinalities of each class in \code{data}, each entry corresponds to one class.
}
  \item{pretransform}{
The method of data scaling. 

\code{NULL} to use the original data, 

The data may be scaled jointly or separately:

\code{1Mom} or \code{1MCD} for joint scaling of the classes, 

\code{NMom} or \code{NMCD} for separate scaling of the classes.

You may use traditional moments or Minimum Covariance Determinant (MCD) estimates for mean and covariance:

\code{1Mom} or \code{NMom} for scaling using traditional data moments, 

\code{1MCD} or \code{NMCD} for scaling using robust MCD data moments.
}
  \item{kernel}{
\code{"EDKernel"} for the kernel of type 1/(1+kernel.bandwidth*EuclidianDistance2(x, y)), 

\code{"GKernel"} [default and recommended] for the simple Gaussian kernel, 

\code{"EKernel"} exponential kernel: exp(-kernel.bandwidth*EuclidianDistance(x, y)), 

%\code{"TriangleKernel"}, 
\code{"VarGKernel"} variable Gaussian kernel, where \code{kernel.bandwidth} is proportional to the \code{depth.zonoid} of a point.
}
  \item{kernel.bandwidth}{
  the bandwidth parameter of the kernel. If \code{NULL} - the Scott's rule of thumb is used.
  May be a single value for all classes, or a vector of values for each of the classes.
}
  \item{mah.parMcd}{
is the value of the argument \code{alpha} for the function \code{\link{covMcd}}; is used when \code{pretransform = "*MCD"}.
}
}
\details{
The potential representation is calculated in the same way as in \code{\link{depth.potential}}, see References below for more information and details.
}
\value{
Matrix of objects, each object (row) is represented via its potentials (columns) w.r.t. each of the classes of the training sample; order of the classes in columns corresponds to the one in the argument \code{cardinalities}.
}
\references{
Aizerman, M.A., Braverman, E.M., and Rozonoer, L.I. (1970). \emph{The Method of Potential Functions in the Theory of Machine Learning}. Nauka (Moscow).

Pokotylo, O. and Mosler, K. (2015). Classification with the pot-pot plot. \emph{Mimeo}.
}
\seealso{
\code{\link{ddalpha.train}} and \code{\link{ddalpha.classify}} for application, \code{\link{depth.potential}} for calculation of the potential.
}
\examples{

# Generate a bivariate normal location-shift classification task
# containing 20 training objects
class1 <- mvrnorm(50, c(0,0), 
                  matrix(c(1,1,1,4), nrow = 2, ncol = 2, byrow = TRUE))
class2 <- mvrnorm(50, c(1,1), 
                  matrix(c(1,1,1,4), nrow = 2, ncol = 2, byrow = TRUE))
data <- rbind(class1, class2)
plot(data, col = c(rep(1,50), rep(2,50)))
# potential with rule of thumb bandwidth
ds = depth.space.potential(data, c(50, 50))
# draw.ddplot(depth.space = ds, cardinalities = c(50, 50))

# potential with bandwidth = 0.5 and joint scaling
ds = depth.space.potential(data, c(50, 50), kernel.bandwidth = 0.5,
                           pretransform = "1Mom")
# draw.ddplot(depth.space = ds, cardinalities = c(50, 50))

# potential with bandwidth = 0.5 and separate scaling
ds = depth.space.potential(data, c(50, 50), kernel.bandwidth = 0.5, 
                           pretransform = "NahMom") # or without pretransform
# draw.ddplot(depth.space = ds, cardinalities = c(50, 50))

data <- getdata("hemophilia")
cardinalities = c(sum(data$gr == "normal"), sum(data$gr == "carrier"))
ds = depth.space.potential(data[,1:2], cardinalities)
# draw.ddplot(depth.space = ds, cardinalities = cardinalities)

}
\keyword{ robust }
\keyword{ multivariate }
\keyword{ nonparametric }
