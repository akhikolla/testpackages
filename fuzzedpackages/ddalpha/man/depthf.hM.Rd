\name{depthf.hM}
\alias{depthf.hM}
\title{h-Mode Depth for Functional Data}
\usage{
depthf.hM(datafA, datafB, range = NULL, d = 101, norm = c("C", "L2"),
  q = 0.2)
}
\arguments{
\item{datafA}{Functions whose depth is computed, represented by a \code{dataf} object of their arguments
and functional values. \code{m} stands for the number of functions.}

\item{datafB}{Random sample functions with respect to which the depth of \code{datafA} is computed. 
\code{datafB} is represented by a \code{dataf} object of their arguments
 and functional values. \code{n} is the sample size. The grid of observation points for the 
functions \code{datafA} and \code{datafB} may not be the same.}

\item{range}{The common range of the domain where the functions \code{datafA} and \code{datafB} are observed.
Vector of length 2 with the left and the right end of the interval. Must contain all arguments given in 
\code{datafA} and \code{datafB}.}

\item{d}{Grid size to which all the functional data are transformed. For depth computation, 
all functional observations are first transformed into vectors of their functional values of length \code{d}
corresponding to equi-spaced points in the domain given by the interval \code{range}. Functional values in these
points are reconstructed using linear interpolation, and extrapolation.}

\item{norm}{The norm used for the computation of the depth. Two possible 
choices are implemented: \code{C} for the uniform norm of continuous functions, 
and \code{L2} for the \eqn{L^2} norm of integrable functions.}

\item{q}{The quantile used to determine the value of the bandwidth \eqn{h}
in the computation of the h-mode depth. \eqn{h} is taken as the \code{q}-quantile of
all non-zero distances between the functions \code{B}. By default, this value is set
to \code{q=0.2}, in accordance with the choice of Cuevas et al. (2007).}
}
\value{
A vector of length \code{m} of the h-mode depth values.
}
\description{
The h-mode depth of functional real-valued data.
}
\details{
The function returns the vectors of the sample h-mode depth values. The kernel 
used in the evaluation is the standard Gaussian kernel, the bandwidth value is chosen
as a quantile of the non-zero distances between the random sample curves.
}
\examples{
datafA = dataf.population()$dataf[1:20]
datafB = dataf.population()$dataf[21:50]
depthf.hM(datafA,datafB)
depthf.hM(datafA,datafB,norm="L2")
}
\references{
Cuevas, A., Febrero, M. and Fraiman, R.  (2007).
Robust estimation and classification for functional data via projection-based depth notions. 
\emph{Computational Statistics} \bold{22} (3), 481--496.

Nagy, S., Gijbels, I. and Hlubinka, D. (2016). 
Weak convergence of discretely observed functional data with applications. 
\emph{Journal of Multivariate Analysis}, \bold{146}, 46--62.
}
\author{
Stanislav Nagy, \email{nagy at karlin.mff.cuni.cz}
}
\keyword{depth}
\keyword{functional}
