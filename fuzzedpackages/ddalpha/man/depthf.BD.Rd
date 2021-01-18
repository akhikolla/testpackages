\name{depthf.BD}
\alias{depthf.BD}
\title{Band Depth for Functional Data}
\usage{
depthf.BD(datafA, datafB, range = NULL, d = 101)
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
}
\value{
A vector of length \code{m} of the band depth values.
}
\description{
The (unadjusted) band depth 
for functional real-valued data of order \code{J=2}.
}
\details{
The function returns the vector of the sample (unadjusted) band depth values.
}
\examples{
datafA = dataf.population()$dataf[1:20]
datafB = dataf.population()$dataf[21:50]
depthf.BD(datafA,datafB)

}
\references{
Lopez-Pintado, S. and Romo, J. (2009), On the concept of depth for functional data,
\emph{J. Amer. Statist. Assoc.} \bold{104} (486), 718 - 734.
}
\seealso{
\code{\link{depthf.ABD}}, \code{\link{depthf.fd1}}
}
\author{
Stanislav Nagy, \email{nagy at karlin.mff.cuni.cz}
}
\keyword{depth}
\keyword{functional}
