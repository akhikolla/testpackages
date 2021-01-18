\name{depthf.fd1}
\alias{depthf.fd1}
\title{Univariate Integrated and Infimal Depth for Functional Data}
\usage{
depthf.fd1(datafA, datafB, range = NULL, d = 101, order = 1, approx = 0)
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

\item{order}{The order of the order extended integrated and infimal depths.
By default, this is set to \code{1}, meaning that the usual univariate depths of 
the functional values are computed. For \code{order=2} or \code{3}, the second
and the third order extended integrated and infimal depths are computed, 
respectively.}

\item{approx}{Number of approximations used in the computation of the order extended depth
for \code{order} greater than \code{1}. For \code{order=2}, the default
value is set to \code{0}, meaning that the depth is computed at all possible \code{d^order}
combinations of the points in the domain. For \code{order=3}, 
the default value is set to \code{101}. When \code{approx} is a positive integer, \code{approx}
points are randomly sampled in \code{[0,1]^order} and at these points the \code{order}-variate depths of the
corresponding functional values are computed.}
}
\value{
Four vectors of length \code{m} of depth values are returned:
\itemize{
\item \code{Simpl_FD} the integrated depth based on the simplicial depth,
\item \code{Half_FD} the integrated depth based on the halfspace depth,
\item \code{Simpl_ID} the infimal depth based on the simplicial depth,
\item \code{Half_ID} the infimal depth based on the halfspace depth.
}
In addition, two vectors of length \code{m} of the relative area of smallest depth values is returned:
\itemize{
\item \code{Simpl_IA} the proportions of points at which the depth \code{Simpl_ID} was attained,
\item \code{Half_IA} the proportions of points at which the depth \code{Half_ID} was attained.
}
The values \code{Simpl_IA} and \code{Half_IA} are always in the interval [0,1]. 
They introduce ranking also among functions having the same
infimal depth value - if two functions have the same infimal depth, the one with larger infimal area
\code{IA} is said to be less central. 	
For \code{order=2} and \code{m=1}, two additional matrices of pointwise depths are also returned:
\itemize{
   \item \code{PSD} the matrix of size \code{d*d} containing the computed 
   pointwise bivariate simplicial depths used for the computation of \code{Simpl_FD} and \code{Simpl_ID},
   \item \code{PHD} the matrix of size \code{d*d} containing the computed 
   pointwise bivariate halfspace depths used for the computation of \code{Half_FD} and \code{Half_ID}.
}
For \code{order=3}, only \code{Half_FD} and \code{Half_ID} are provided.
}
\description{
Usual, and order extended integrated and infimal depths for real-valued functional data based on the
halfspace and simplicial depth.
}
\details{
The function returns vectors of sample integrated and infimal depth values.
}
\examples{
datafA = dataf.population()$dataf[1:20]
datafB = dataf.population()$dataf[21:50]
depthf.fd1(datafA,datafB)
depthf.fd1(datafA,datafB,order=2)
depthf.fd1(datafA,datafB,order=3,approx=51)

}
\references{
Nagy, S., Gijbels, I. and Hlubinka, D. (2016). 
Weak convergence of discretely observed functional data with applications. 
\emph{Journal of Multivariate Analysis}, \bold{146}, 46--62.

Nagy, S., Gijbels, I. and Hlubinka, D.  (2017).
Depth-based recognition of shape outlying functions. 
\emph{Journal of Computational and Graphical Statistics}, \bold{26} (4), 883--893.
}
\seealso{
\code{\link{depthf.fd2}}, \code{\link{infimalRank}}
}
\author{
Stanislav Nagy, \email{nagy at karlin.mff.cuni.cz}
}
\keyword{depth}
\keyword{functional}
