\name{depthf.ABD}
\alias{depthf.ABD}
\title{Adjusted Band Depth for Functional Data}
\usage{
depthf.ABD(datafA, datafB, range = NULL, d = 101, norm = c("C", "L2"),
  J = 2, K = 1)
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
points are reconstructed using linear interpolation, and extrapolation, see Nagy et al. (2016).}

\item{norm}{The norm used for the computation of the depth. Two possible 
choices are implemented: \code{C} for the uniform norm of continuous functions, 
and \code{L2} for the \eqn{L^2} norm of integrable functions.}

\item{J}{The order of the adjusted band depth, that is the maximal number of functions
taken in a band. Acceptable values are \code{2}, \code{3},... By default this value is set to \code{2}. 
Note that this is NOT the order as
defined in the order-extended version of adjusted band depths in Nagy et al. (2016), used
for the detection of shape outlying curves.}

\item{K}{Number of sub-samples of the functions from \code{B} taken to speed up the
computation. By default, sub-sampling is not performed. Values of \code{K} larger than \code{1}
result in an approximation of the adjusted band depth.}
}
\value{
A vectors of length \code{m} of the adjusted band depths.
}
\description{
The adjusted band depth 
of functional real-valued data based on either the
\eqn{C} (uniform) norm, or on the \eqn{L^2} norm of functions.
}
\details{
The function returns the vector of the sample adjusted band depth values. The kernel 
used in the evaluation is the function \eqn{K(u) = exp(-u)}.
}
\examples{
datafA = dataf.population()$dataf[1:20]
datafB = dataf.population()$dataf[21:50]
depthf.ABD(datafA,datafB)
depthf.ABD(datafA,datafB,norm="L2")

}
\references{
Gijbels, I., Nagy, S. (2015).
Consistency of non-integrated depths for functional data.
\emph{Journal of Multivariate Analysis} \bold{140}, 259--282.

Nagy, S., Gijbels, I. and Hlubinka, D. (2016). 
Weak convergence of discretely observed functional data with applications. 
\emph{Journal of Multivariate Analysis}, \bold{146}, 46--62.

Nagy, S., Gijbels, I. and Hlubinka, D.  (2017).
Depth-based recognition of shape outlying functions. 
\emph{Journal of Computational and Graphical Statistics}, \bold{26} (4), 883--893.
}
\seealso{
\code{\link{depthf.BD}}
}
\author{
Stanislav Nagy, \email{nagy at karlin.mff.cuni.cz}
}
\keyword{depth}
\keyword{functional}
