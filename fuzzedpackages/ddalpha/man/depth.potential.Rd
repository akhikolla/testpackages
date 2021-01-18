\name{depth.potential}
\alias{depth.potential}
\title{
Calculate Potential of the Data
}
\description{
Calculate the potential of the points w.r.t. a multivariate data set. The potential is the kernel-estimated density multiplied by the prior probability of a class. Different from the data depths, a density estimate measures at a given point how much mass is located around it.
}
\usage{
depth.potential (x, data, pretransform = "1Mom", 
                kernel = "GKernel", kernel.bandwidth = NULL, mah.parMcd = 0.75)
}
\arguments{
  \item{x}{
Matrix of objects (numerical vector as one object) whose depth is to be calculated; each row contains a \eqn{d}-variate point. Should have the same dimension as \code{data}.
}
  \item{data}{
Matrix of data where each row contains a \eqn{d}-variate point, w.r.t. which the depth is to be calculated.
}
  \item{pretransform}{
The method of data scaling. 

\code{NULL} to use the original data, 

\code{1Mom} or \code{NMom} for scaling using data moments, 

\code{1MCD} or \code{NMCD} for scaling using robust data moments (Minimum Covariance Determinant (MCD) ).
}
  \item{kernel}{
\code{"EDKernel"} for the kernel of type 1/(1+kernel.bandwidth*EuclidianDistance2(x, y)), 

\code{"GKernel"} [default and recommended] for the simple Gaussian kernel, 

\code{"EKernel"} exponential kernel: exp(-kernel.bandwidth*EuclidianDistance(x, y)), 

%\code{"TriangleKernel"}, 
\code{"VarGKernel"} variable Gaussian kernel, where \code{kernel.bandwidth} is proportional to the \code{depth.zonoid} of a point.
}
  \item{kernel.bandwidth}{
  the single bandwidth parameter of the kernel. If \code{NULL} - the Scott's rule of thumb is used.
}
  \item{mah.parMcd}{
is the value of the argument \code{alpha} for the function \code{\link{covMcd}}; is used when \code{pretransform = "*MCD"}.
}
}
\details{
The potential is the kernel-estimated density multiplied by the prior probability of a class.
The kernel bandwidth matrix is decomposed into two parts, one of which describes the form of the data, and the other the width of the kernel. Then the first part is used to transform the data using the moments, while the second is employed as a parameter of the kernel and tuned to achieve the best separation.
For details see Pokotylo and Mosler (2015).
}
\value{
Numerical vector of potentials, one for each row in \code{x}; or one potential value if \code{x} is a numerical vector.
}
\references{
Aizerman, M.A., Braverman, E.M., and Rozonoer, L.I. (1970). \emph{The Method of Potential Functions in the Theory of Machine Learning}. Nauka (Moscow).

Pokotylo, O. and Mosler, K. (2015). Classification with the pot-pot plot. \emph{Mimeo}.
}
\seealso{
\code{\link{depth.halfspace}} for calculation of the Tukey depth.

\code{\link{depth.Mahalanobis}} for calculation of Mahalanobis depth.

\code{\link{depth.projection}} for calculation of projection depth.

\code{\link{depth.simplicial}} for calculation of simplicial depth.

\code{\link{depth.simplicialVolume}} for calculation of simplicial volume depth.

\code{\link{depth.spatial}} for calculation of spatial depth.

\code{\link{depth.zonoid}} for calculation of zonoid depth.
}
\examples{
# 3-dimensional normal distribution
data <- mvrnorm(200, rep(0, 3), 
                matrix(c(1, 0, 0,
                         0, 2, 0, 
                         0, 0, 1),
                       nrow = 3))
x <- mvrnorm(10, rep(1, 3), 
             matrix(c(1, 0, 0,
                      0, 1, 0, 
                      0, 0, 1),
                    nrow = 3))

# potential with rule of thumb bandwidth
pot <- depth.potential(x, data)
cat("Potentials: ", pot, "\n")

# potential with bandwidth = 0.1
pot <- depth.potential(x, data, kernel.bandwidth = 0.1)
cat("Potentials: ", pot, "\n")

# potential with robust MCD scaling
pot <- depth.potential(x, data, kernel.bandwidth = 0.1, 
                      pretransform = "NMCD", mah.parMcd = 0.6)
cat("Potentials: ", pot, "\n")
}
\keyword{ robust }
\keyword{ multivariate }
\keyword{ nonparametric }
