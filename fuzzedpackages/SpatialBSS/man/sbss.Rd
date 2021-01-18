\name{sbss}
\alias{sbss}
\alias{sbss.default}
\alias{sbss.SpatialPointsDataFrame}
\alias{sbss.sf}

\title{
Spatial Blind Source Separation
}

\description{
\code{sbss} estimates the unmixing matrix assuming a spatial blind source separation model by simultaneous/jointly diagonalizing the covariance matrix and one/many local covariance matrices. These local covariance matrices are determined by spatial kernel functions. Three types of such kernel functions are supported.
}

\usage{
sbss(x, ...)

\method{sbss}{default}(x, coords, kernel_type = c('ring', 'ball', 'gauss'), 
     kernel_parameters, lcov = c('lcov', 'ldiff'), 
     ordered = TRUE, kernel_list = NULL, \dots)
\method{sbss}{SpatialPointsDataFrame}(x, \dots)
\method{sbss}{sf}(x, \dots)
}

\arguments{
  \item{x}{
either a numeric matrix of dimension \code{c(n, p)} where the \code{p} columns correspond to the entries of the random field and the \code{n} rows are the observations, an object of class \code{\link[sp]{SpatialPointsDataFrame}} or an object of class \code{\link[sf]{sf}}. 
}
  \item{coords}{
a numeric matrix of dimension \code{c(n,2)} where each row represents the coordinates of a point in the spatial domain. Only needed if \code{x} is a matrix and the argument \code{kernel_list} is \code{NULL}.
}
  \item{kernel_type}{
a string indicating which kernel function to use. Either \code{'ring'} (default), \code{'ball'} or \code{'gauss'}.
}
  \item{kernel_parameters}{
a numeric vector that gives the parameters for the kernel function. At least length of one for \code{'ball'} and \code{'gauss'} or two for \code{'ring'} kernel, see details.
}
  \item{lcov}{
a string indicating which type of local covariance matrix to use. Either \code{'lcov'} (default) or \code{'ldiff'}.
}
  \item{ordered}{
logical. If \code{TRUE} the entries of the latent field are ordered by the sum of squared (pseudo-)eigenvalues of the diagonalized local covariance matrix/matrices. Default is \code{TRUE}.
}
  \item{kernel_list}{
a list of spatial kernel matrices with dimension \code{c(n,n)}, see details. Usually computed by the function \code{\link{spatial_kernel_matrix}}.
}
  \item{\dots}{
further arguments for the fast real joint diagonalization algorithm that jointly diagonalizes the local covariance matrices. See details and \code{\link[JADE]{frjd}}.
}
}

\details{
Two versions of local covariance matrices are implemented, the argument \code{lcov} determines which version is used:
\itemize{
  \item \code{'lcov'}: \deqn{ LCov(f) = 1/n \sum_{i,j} f(d_{i,j}) x(s_i) x(s_j)' ,}
  \item \code{'ldiff'}: \deqn{ LDiff(f) = 1/n \sum_{i,j} f(d_{i,j}) (x(s_i)-x(s_j)) (x(s_i)-x(s_j))' .}
}
Where \eqn{d_{i,j} \ge 0} correspond to the pairwise distances between coordinates, \eqn{x(s_i)} are the \code{p} random field values at location \eqn{s_i} and the kernel function \eqn{f(d)} determines the locality. LDiff matrices are supposed to be more robust when the random field shows a smooth trend. The following kernel functions are implemented and chosen with the argument \code{kernel_type}:
\itemize{
  \item \code{'ring'}: parameters are inner radius \eqn{r_i} and outer radius \eqn{r_o}, with \eqn{r_i < r_o}, and \eqn{r_i, r_o \ge 0}: \deqn{f(d;r_i, r_o) = I(r_i < d \le r_o)}
  \item \code{'ball'}: parameter is the radius \eqn{r}, with \eqn{r \ge 0}: \deqn{f(d;r) = I(d \le r)}
  \item \code{'gauss'}: Gaussian function where 95\% of the mass is inside the parameter \eqn{r}, with \eqn{r \ge 0}: \deqn{f(d;r) = exp(-0.5 (\Phi^{-1}(0.95) d/r)^2)}
}

The argument \code{kernel_type} determines the used kernel function as presented above, the argument \code{kernel_parameters} gives the corresponding parameters for the kernel function. Specifically, if \code{kernel_type} equals \code{'ball'} or \code{'gauss'} then \code{kernel_parameters} is a numeric vector where each entry corresponds to one parameter. Hence, \code{length(kernel_parameters)} local covariance matrices are used. Whereas, if \code{kernel_type} equals \code{'ring'}, then \code{kernel_parameters} must be a numeric vector of even length where subsequently the inner and outer radii must be given (informally: \code{c(r_i1, r_o1, r_i2, r_o2, \dots)}). In that case \code{length(kernel_parameters) / 2} local covariance matrices are used.

Internally, \code{sbss} calls \code{\link{spatial_kernel_matrix}} to compute a list of \code{c(n,n)} kernel matrices based on the parameters given, where each entry of those matrices corresponds to \eqn{f(d_{i,j})}. Alternatively, such a list of kernel matrices can be given directly to the function \code{sbss} via the \code{kernel_list} argument. This is useful when \code{sbss} is called numerous times with the same coordinates/kernel functions as the computation of the kernel matrices is then done only once prior the actual \code{sbss} calls. For details see also \code{\link{spatial_kernel_matrix}}.

If more than one local covariance matrix is used \code{sbss} jointly diagonalizes these matrices with the function \code{\link[JADE]{frjd}}. \code{\dots} provides arguments for \code{frjd}, useful arguments might be:
\itemize{
  \item \code{eps}: tolerance for convergence.
  \item \code{maxiter}: maximum number of iterations.
}

}

\value{
\code{sbss} returns a list of class \code{'sbss'} with the following entries: 
  \item{s}{
  object of \code{class(x)} containing the estimated source random field.
}  
  \item{coords}{
  coordinates of the observations. Is \code{NULL} if \code{x} was a matrix and the argument \code{kernel_list} was not \code{NULL} at the \code{sbss} call.
}

  \item{w}{
  estimated unmixing matrix.
}

  \item{w_inv}{
  inverse of the estimated unmixing matrix.
}

  \item{d}{
  matrix of stacked (jointly) diagonalized local covariance matrices with dimension \code{c(length(kernel_parameters)*p,p)} for \code{'ball'} and \code{'gauss'} kernel or \code{c(}
  \code{(length(kernel_parameters) / 2)*p,p)} for \code{'ring'} kernel.
}

  \item{x_mu}{
  columnmeans of \code{x}.
}

  \item{cov_inv_sqrt}{
  square root of the inverse sample covariance matrix of \code{x}.
}

}

\references{
Bachoc, F., Genton, M. G, Nordhausen, K., Ruiz-Gazen, A. and Virta, J. (2020), \emph{Spatial Blind Source Separation}, Biometrika, \url{https://doi.org/10.1093/biomet/asz079}, \cr <doi:10.1093/biomet/asz079>.

Nordhausen, K., Oja, H., Filzmoser, P., Reimann, C. (2015), \emph{Blind Source Separation for Spatial Compositional Data}, Mathematical Geosciences 47, 753-770, \url{https://doi.org/10.1007/s11004-014-9559-5}, <https://doi.org/10.1007/s11004-014-9559-5>.

}

\seealso{
\code{\link{spatial_kernel_matrix}}, \code{\link{local_covariance_matrix}}, \code{\link[sp]{sp}}, \code{\link[sf]{sf}}, \code{\link[JADE]{frjd}}
}

\examples{
# simulate coordinates
coords <- runif(1000 * 2) * 20
dim(coords) <- c(1000, 2)

# simulate random field
if (!requireNamespace('RandomFields', quietly = TRUE)) {
  stop('Please install the package RandomFields to run the example code.')
} else {
  RandomFields::RFoptions(spConform = FALSE)
  field_1 <- RandomFields::RFsimulate(model = RandomFields::RMexp(), 
                                      x = coords)
  field_2 <- RandomFields::RFsimulate(model = RandomFields::RMspheric(), 
                                      x = coords)
  field_3 <- RandomFields::RFsimulate(model = RandomFields::RMwhittle(nu = 2),  
                                      x = coords)
  field <- cbind(field_1, field_2, field_3)
}

# apply sbss with three ring kernels
kernel_parameters <- c(0, 1, 1, 2, 2, 3)
sbss_result <- 
  sbss(field, coords, kernel_type = 'ring', kernel_parameters = kernel_parameters)

# print object
print(sbss_result)

# plot latent field
plot(sbss_result, colorkey = TRUE, as.table = TRUE, cex = 1)

# predict latent fields on grid
predict(sbss_result, colorkey = TRUE, as.table = TRUE, cex = 1)

# unmixing matrix
w_unmix <- coef(sbss_result)

# apply the same sbss with a kernel list
kernel_list <- spatial_kernel_matrix(coords, kernel_type = 'ring', kernel_parameters)
sbss_result_k <- sbss(field, kernel_list = kernel_list)

# apply sbss with three ring kernels and local difference matrices
sbss_result_ldiff <- 
  sbss(field, coords, kernel_type = 'ring', 
       kernel_parameters = kernel_parameters, lcov = 'ldiff')


}

\keyword{ multivariate }
\keyword{ spatial }
