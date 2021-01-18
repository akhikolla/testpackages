\name{spatial_kernel_matrix}
\alias{spatial_kernel_matrix}

\title{
Computation of Spatial Kernel Matrices
}

\description{
\code{spatial_kernel_matrix} computes spatial kernel matrices for a given kernel function with its parameters and a set of coordinates.
}

\usage{
spatial_kernel_matrix(coords, kernel_type = c('ring', 'ball', 'gauss'), 
                      kernel_parameters)
}

\arguments{
  \item{coords}{
a numeric matrix of dimension \code{c(n,2)} where each row represents the coordinates of a point in the spatial domain.
}
  \item{kernel_type}{
a character string indicating which kernel function to use. Either \code{'ring'} (default), \code{'ball'} or \code{'gauss'}.
}
  \item{kernel_parameters}{
a numeric vector that gives the parameters for the kernel function. At least length of one for \code{'ball'} and \code{'gauss'} or two for \code{'ring'} kernel, see details.
}
}

\details{
Two versions of local covariance matrices can be defined:
\itemize{
  \item \code{'lcov'}: \deqn{ LCov(f) = 1/n \sum_{i,j} f(d_{i,j}) x(s_i) x(s_j)' ,}
  \item \code{'ldiff'}: \deqn{ LDiff(f) = 1/n \sum_{i,j} f(d_{i,j}) (x(s_i)-x(s_j)) (x(s_i)-x(s_j))' .}
}
Where \eqn{d_{i,j} \ge 0} correspond to the pairwise distances between coordinates, \eqn{x(s_i)} are the \code{p} random field values at location \eqn{s_i} and the kernel function \eqn{f(d)} determines the locality. The function \code{spatial_kernel_matrix} computes a list of \code{c(n,n)} matrices where each entry of these matrices correspond to the spatial kernel function evaluated at the distance between two points, mathematically the entry ij of each kernel matrix is \eqn{f(d_{i,j})}. The following kernel functions are implemented and chosen with the argument \code{kernel_type}:
\itemize{
  \item \code{'ring'}: parameters are inner radius \eqn{r_i} and outer radius \eqn{r_o}, with \eqn{r_i < r_o}, and \eqn{r_i, r_o \ge 0}: \deqn{f(d;r_i, r_o) = I(r_i < d \le r_o)}
  \item \code{'ball'}: parameter is the radius \eqn{r}, with \eqn{r \ge 0}: \deqn{f(d;r) = I(d \le r)}
  \item \code{'gauss'}: Gaussian function where 95\% of the mass is inside the parameter \eqn{r}, with \eqn{r \ge 0}: \deqn{f(d;r) = exp(-0.5 (\Phi^{-1}(0.95) d/r)^2)}
}

The argument \code{kernel_type} determines the used kernel function as presented above, the argument \code{kernel_parameters} gives the corresponding parameters for the kernel function. Specifically, if \code{kernel_type} equals \code{'ball'} or \code{'gauss'} then \code{kernel_parameters} is a numeric vector where each entry corresponds to one parameter. Hence, \code{length(kernel_parameters)} spatial kernel matrices of type \code{kernel_type} are computed. Whereas, if \code{kernel_type} equals \code{'ring'}, then \code{kernel_parameters} must be a numeric vector of even length where subsequently the inner and outer radii must be given (informally: \code{c(r_i1, r_o1, r_i2, r_o2, \dots)}). In that case \code{length(kernel_parameters) / 2} spatial kernel matrices of type \code{'ring'} are computed.

The output of this function can be used with the function \code{\link{sbss}} to avoid unnecessary computation of kernel matrices when \code{\link{sbss}} is called multiple times with the same coordinate/kernel function setting. Additionally, the output can be used with the function \code{\link{local_covariance_matrix}} to actually compute local covariance matrices as defined above based on a given set of spatial kernel matrices.
}

\value{
\code{spatial_kernel_matrix} returns a list with length of \code{length(kernel_parameters)} (for \code{'ball'} and \code{'gauss'} kernel functions) or \code{length(kernel_parameters) / 2} (for \code{'ring'} kernel function) containing numeric matrices of dimension \code{c(n,n)} corresponding to the spatial kernel matrices.
}

\references{
Bachoc, F., Genton, M. G, Nordhausen, K., Ruiz-Gazen, A. and Virta, J. (2020), \emph{Spatial Blind Source Separation}, Biometrika, \url{https://doi.org/10.1093/biomet/asz079}, \cr <doi:10.1093/biomet/asz079>.
}

\seealso{
\code{\link{sbss}}, \code{\link{local_covariance_matrix}}
}

\examples{
# simulate a set of coordinates
coords <- rnorm(100 * 2)
dim(coords) <- c(100, 2)

# computing two ring kernel matrices
kernel_params_ring <- c(0, 0.5, 0.5, 2)
ring_kernel_list <- 
  spatial_kernel_matrix(coords, 'ring', kernel_params_ring)

# computing three ball kernel matrices
kernel_params_ball <- c(0.5, 1, 2)
ball_kernel_list <- 
  spatial_kernel_matrix(coords, 'ball', kernel_params_ball)

# computing three gauss kernel matrices
kernel_params_gauss <- c(0.5, 1, 2)
gauss_kernel_list <- 
  spatial_kernel_matrix(coords, 'gauss', kernel_params_gauss)

}

\keyword{ array }

