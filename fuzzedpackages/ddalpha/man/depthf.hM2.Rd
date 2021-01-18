\name{depthf.hM2}
\alias{depthf.hM2}
\title{Bivariate h-Mode Depth for Functional Data Based on the \eqn{L^2} Metric}
\usage{
depthf.hM2(datafA, datafB, range = NULL, d = 101, q = 0.2)
}
\arguments{
\item{datafA}{Bivariate functions whose depth is computed, represented by a multivariate \code{dataf} object of 
their arguments (vector), and a matrix with two columns of the corresponding bivariate functional values. 
\code{m} stands for the number of functions.}

\item{datafB}{Bivariate random sample functions with respect to which the depth of \code{datafA} is computed. 
\code{datafB} is represented by a multivariate \code{dataf} object of their arguments
 (vector), and a matrix with two columns of the corresponding bivariate functional values.
\code{n} is the sample size. The grid of observation points for the 
functions \code{datafA} and \code{datafB} may not be the same.}

\item{range}{The common range of the domain where the functions \code{datafA} and \code{datafB} are observed.
Vector of length 2 with the left and the right end of the interval. Must contain all arguments given in 
\code{datafA} and \code{datafB}.}

\item{d}{Grid size to which all the functional data are transformed. For depth computation, 
all functional observations are first transformed into vectors of their functional values of length \code{d}
corresponding to equi-spaced points in the domain given by the interval \code{range}. Functional values in these
points are reconstructed using linear interpolation, and extrapolation.}

\item{q}{The quantile used to determine the value of the bandwidth \eqn{h}
in the computation of the h-mode depth. \eqn{h} is taken as the \code{q}-quantile of
all non-zero distances between the functions \code{B}. By default, this value is set
to \code{q=0.2}, in accordance with the choice of Cuevas et al. (2007).}
}
\value{
Three vectors of length \code{m} of h-mode depth values are returned:
\itemize{
\item \code{hM} the unscaled h-mode depth,
\item \code{hM_norm} the h-mode depth \code{hM} linearly transformed so that its range is [0,1],
\item \code{hM_norm2} the h-mode depth \code{FD} linearly transformed by a transformation such that 
   the range of the h-mode depth of \code{B} with respect to \code{B} is [0,1]. This depth may give negative values.
}
}
\description{
The h-mode depth 
of functional bivariate data (that is, data of the form \eqn{X:[a,b] \to R^2},
or \eqn{X:[a,b] \to R} and the derivative of \eqn{X}) based on the
\eqn{L^2} metric of functions.
}
\details{
The function returns the vectors of sample h-mode depth values. The kernel 
used in the evaluation is the standard Gaussian kernel, the bandwidth value is chosen
as a quantile of the non-zero distances between the random sample curves.
}
\examples{
datafA = dataf.population()$dataf[1:20]
datafB = dataf.population()$dataf[21:50]

datafA2 = derivatives.est(datafA,deriv=c(0,1))
datafB2 = derivatives.est(datafB,deriv=c(0,1))

depthf.hM2(datafA2,datafB2)

depthf.hM2(datafA2,datafB2)$hM
# depthf.hM2(cbind(A2[,,1],A2[,,2]),cbind(B2[,,1],B2[,,2]))$hM
# the two expressions above should give the same result

}
\references{
Cuevas, A., Febrero, M. and Fraiman, R.  (2007). 
Robust estimation and classification for functional data via projection-based depth notions. 
\emph{Computational Statistics} \bold{22} (3), 481--496.
}
\seealso{
\code{\link{depthf.hM}}
}
\author{
Stanislav Nagy, \email{nagy at karlin.mff.cuni.cz}
}
\keyword{depth}
\keyword{derivatives}
\keyword{functional}
