\name{depthf.RP1}
\alias{depthf.RP1}
\title{Univariate Random Projection Depths for Functional Data}
\usage{
depthf.RP1(datafA, datafB, range = NULL, d = 101, nproj = 50, nproj2 = 5)
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

\item{nproj}{Number of projections taken in the computation of the random projection depth. By default taken
to be \code{51}.}

\item{nproj2}{Number of projections taken in the computation of the random functional depth. By default taken
to be \code{5}. \code{nproj2} should be much smaller than \code{d}, the dimensionality of the discretized 
functional data.}
}
\value{
Three vectors of depth values of length \code{m} are returned:
\itemize{
\item \code{Simpl_FD} the random projection depth based on the univariate simplicial depth,
\item \code{Half_FD} the random projection depth based on the univariate halfspace depth,
\item \code{RHalf_FD} the random halfspace depth.
}
}
\description{
Random projection depth and random functional depth for functional data.
}
\details{
The function returns the vectors of sample random projection, and random functional depth values. 
The random projection depth described in Cuevas et al. (2007) is based on the average univariate depth
of one-dimensional projections of functional data. The projections are taken randomly as a sample of standard
normal \code{d}-dimensional random variables, where \code{d} stands for the dimensionality of the discretized
functional data. 

The random functional depth (also called random Tukey depth, or random halfspace depth) is described in
Cuesta-Albertos and Nieto-Reyes (2008). The functional data are projected into the real line in random 
directions as for the random projection depths. Afterwards, an approximation of the halfspace (Tukey) depth
based on this limited number of univariate projections is assessed.
}
\examples{
datafA = dataf.population()$dataf[1:20]
datafB = dataf.population()$dataf[21:50]

depthf.RP1(datafA,datafB)

}
\references{
Cuevas, A., Febrero, M. and Fraiman, R.  (2007).
Robust estimation and classification for functional data via projection-based depth notions, 
\emph{Computational Statistics} \bold{22} (3), 481--496.

Cuesta-Albertos, J.A. and Nieto-Reyes, A. (2008).
 The random Tukey depth. 
\emph{Computational Statistics & Data Analysis} \bold{52} (11), 4979--4988.
}
\seealso{
\code{\link{depthf.RP2}}
}
\author{
Stanislav Nagy, \email{nagy at karlin.mff.cuni.cz}
}
\keyword{depth}
\keyword{functional}
