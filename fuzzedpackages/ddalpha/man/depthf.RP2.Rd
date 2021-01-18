\name{depthf.RP2}
\alias{depthf.RP2}
\title{Bivariate Random Projection Depths for Functional Data}
\usage{
depthf.RP2(datafA, datafB, range = NULL, d = 101, nproj = 51)
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

\item{nproj}{Number of projections taken in the computation of the double random projection depth. By default taken
to be \code{51}.}
}
\value{
Five vectors of length \code{m} are returned:
\itemize{
\item \code{Simpl_FD} the double random projection depth RP2 based on the bivariate simplicial depth,
\item \code{Half_FD} the double random projection depth RP2 based on the bivariate halfspace depth,
\item \code{hM_FD} the double random projection depth RP2 based on the bivariate h-mode depth,
\item \code{Simpl_DD} the double random projection depth RPD based on the univariate simplicial depth,
\item \code{Half_DD} the random projection depth RPD based on the univariate halfspace depth,
}
}
\description{
Double random projection depths of functional bivariate data (that is, data of the form \eqn{X:[a,b] \to R^2},
or \eqn{X:[a,b] \to R} and the derivative of \eqn{X}).
}
\details{
The function returns the vectors of sample double random projection depth values. 
The double random projection depths are described in Cuevas et al. (2007). They are of two types: RP2 type, and
RPD type. Both types of depths are based on bivariate projections of the bivariate functional data. 
These projections are taken randomly as a sample of standard
normal \code{d}-dimensional random variables, where \code{d} stands for the dimensionality of the internally 
represented discretized
functional data. For RP2 type depths, the average bivariate depth of the projected quantities is assessed.
For RPD type depths, further univariate projections of these bivariate projected quantities are evaluated, and
based on these final univariate quantities, the average univariate depth is computed.
}
\examples{
datafA = dataf.population()$dataf[1:20]
datafB = dataf.population()$dataf[21:50]

dataf2A = derivatives.est(datafA,deriv=c(0,1))
dataf2B = derivatives.est(datafB,deriv=c(0,1))
depthf.RP2(dataf2A,dataf2B)


}
\references{
Cuevas, A., Febrero, M. and Fraiman, R.  (2007).
Robust estimation and classification for functional data via projection-based depth notions.
\emph{Computational Statistics} \bold{22} (3), 481--496.
}
\seealso{
\code{\link{depthf.RP1}}
}
\author{
Stanislav Nagy, \email{nagy at karlin.mff.cuni.cz}
}
\keyword{depth}
\keyword{derivatives}
\keyword{functional}
