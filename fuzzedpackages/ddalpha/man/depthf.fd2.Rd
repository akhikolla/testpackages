\name{depthf.fd2}
\alias{depthf.fd2}
\title{Bivariate Integrated and Infimal Depth for Functional Data}
\usage{
depthf.fd2(datafA, datafB, range = NULL, d = 101)
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
}
\value{
Four vectors of length \code{m} are returned:
\itemize{
\item \code{Simpl_FD} the integrated depth based on the bivariate simplicial depth,
\item \code{Half_FD} the integrated depth based on the bivariate halfspace depth,
\item \code{Simpl_ID} the infimal depth based on the bivariate simplicial depth,
\item \code{Half_ID} the infimal depth based on the bivariate halfspace depth.
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
}
\description{
Integrated and infimal depths 
of functional bivariate data (that is, data of the form \eqn{X:[a,b] \to R^2},
or \eqn{X:[a,b] \to R} and the derivative of \eqn{X}) based on the
bivariate halfspace and simplicial depths.
}
\details{
The function returns the vectors of sample integrated and infimal depth values.
}
\examples{
datafA = dataf.population()$dataf[1:20]
datafB = dataf.population()$dataf[21:50]

dataf2A = derivatives.est(datafA,deriv=c(0,1))
dataf2B = derivatives.est(datafB,deriv=c(0,1))
depthf.fd2(dataf2A,dataf2B)

}
\references{
Hlubinka, D., Gijbels, I., Omelka, M. and Nagy, S. (2015). 
Integrated data depth for smooth functions and its application in supervised classification. 
\emph{Computational Statistics}, \bold{30} (4), 1011--1031.

Nagy, S., Gijbels, I. and Hlubinka, D.  (2017).
Depth-based recognition of shape outlying functions. 
\emph{Journal of Computational and Graphical Statistics}, \bold{26} (4), 883--893.
}
\seealso{
\code{\link{depthf.fd1}}, \code{\link{infimalRank}}
}
\author{
Stanislav Nagy, \email{nagy at karlin.mff.cuni.cz}
}
\keyword{depth}
\keyword{derivatives}
\keyword{functional}
