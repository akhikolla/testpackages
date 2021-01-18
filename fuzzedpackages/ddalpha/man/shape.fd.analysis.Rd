\name{shape.fd.analysis}
\alias{shape.fd.analysis}
\title{Diagnostic Plot for First and Second Order Integrated and Infimal Depths}
\usage{
shape.fd.analysis(datafA, datafB, range = NULL, d = 101, order = 1,
  method = c("halfspace", "simplicial"), approx = 0, title = "",
  nfun = 10, plot = TRUE)
}
\arguments{
\item{datafA}{A single function whose depth is computed, represented by a 
\code{dataf} object of arguments and functional values.}

\item{datafB}{Functional dataset with respect to which the depth of \code{datafA} is computed. 
\code{datafB} is represented by a \code{dataf} object of arguments and functional values. 
\code{n} stands for the number of functions. The grid of observation points for the 
functions in \code{datafA} and \code{datafB} may not be the same.}

\item{range}{The common range of the domain where the functions \code{datafA} and \code{datafB} are observed.
Vector of length 2 with the left and the right end of the interval. Must contain all arguments given in 
\code{datafA} and \code{datafB}.}

\item{d}{Grid size to which all the functional data are transformed. For depth computation, 
all functional observations are first transformed into vectors of their functional values of length \code{d}
corresponding to equi-spaced points in the domain given by the interval \code{range}. Functional values in these
points are reconstructed using linear interpolation, and extrapolation.}

\item{order}{The order of the depth to be used in the plot, for \code{order=1} produces
the plot of univariate marginal depth of \code{A} and \code{nfun} functions from \code{B} 
over the domain of the functions. For \code{order=2} produces the bivariate contour plot 
of the bivariate depths of \code{A} at couples of points from the domain.}

\item{method}{The depth that is used in the diagnostic plot. possible values are \code{halfspace} for 
the halfspace depth, or \code{simplicial} for the simplicial depth.}

\item{approx}{For \code{order=2}, the number of approximations used in the computation of the order extended depth. By default
this is set to \code{0}, meaning that the depth is computed at all possible \code{d^2}
combinations of the points in the domain. When set to a positive integer, \code{approx}
bivariate points are randomly sampled in unit square, and at these points the bivariate depths of the
corresponding functional values are computed.}

\item{title}{The title of the diagnostic plot.}

\item{nfun}{For \code{order=1}, the number of functions from \code{B} whose coordinate-wise
univariate depths of functional values should be displayed with the depth of \code{A}.
The depth of \code{A} is displayed in solid red line, the depths of the functions from \code{B}
in dashed black.}

\item{plot}{Logical: should the function by plotted?}
}
\value{
For \code{order=1} two depth values, and two vectors of pointwise depths:
\itemize{
\item \code{Simpl_FD} the first order integrated depth based on the simplicial depth,
\item \code{Half_FD} the first order integrated depth based on the halfspace depth,
\item \code{Simpl_ID} the first order infimal depth based on the simplicial depth,
\item \code{Half_ID} the first order infimal depth based on the halfspace depth,
   \item \code{PSD} the vector of length \code{d} containing the computed 
   pointwise univariate simplicial depths used for the computation of \code{Simpl_FD} and \code{Simpl_ID},
   \item \code{PHD} the vector of length \code{d} containing the computed 
   pointwise univariate halfspace depths used for the computation of \code{Half_FD} and \code{Half_ID}.
}
   In addition, the first order integrated / infimal depth diagnostic plot of the function \code{A} with respect to
   the random sample given by the functions corresponding to the rows of the matrix \code{B} is produced.

   For \code{order=2} four depth values, and two matrices of pointwise depths:
\itemize{
\item \code{Simpl_FD} the second order integrated depth based on the simplicial depth,
\item \code{Half_FD} the second order integrated depth based on the halfspace depth,
\item \code{Simpl_ID} the second order infimal depth based on the simplicial depth,
\item \code{Half_ID} the second order infimal depth based on the halfspace depth,
   \item \code{PSD} the matrix of size \code{d*d} containing the computed 
   pointwise bivariate simplicial depths used for the computation of \code{Simpl_FD} and \code{Simpl_ID},
   \item \code{PHD} the matrix of size \code{d*d} containing the computed 
   pointwise bivariate halfspace depths used for the computation of \code{Half_FD} and \code{Half_ID}.
}
   In addition, the second order integrated / infimal depth diagnostic plot of the function \code{A} with respect to
   the random sample given by the functions corresponding to the rows of the matrix \code{B} is produced.
}
\description{
Produce the diagnostic plot based on the fist or second order extended integrated / infimal depths.
}
\details{
Plots a diagnostic plot of pointwise univariate (or bivariate) depths for all possible points (or couples of points) from the domain of the 
functional data. From such a plot it is possible to infer into the first order (or second order) properties of a single function \emph{x} with respect 
to the given set of functional data. For \code{order=1}, the integral of the displayed function is the integrated depth of \emph{x}, 
the smallest value of the function is the infimal depth of \emph{x}. 
For \code{order=2}, the bivariate integral of the displayed surface gives the second order extended 
integrated depth of \emph{x}, the infimum of this bivariate function gives the second order infimal depth of \emph{x}. 
For details see Nagy et al. (2016) and \code{\link{depthf.fd1}}.
}
\examples{
datafA = dataf.population()$dataf[1]
dataf = dataf.population()$dataf[2:20]
shape.fd.analysis(datafA,dataf,order=1)
shape.fd.analysis(datafA,dataf,order=2,approx=0)

}
\references{
Nagy, S., Gijbels, I. and Hlubinka, D.  (2017).
Depth-based recognition of shape outlying functions. 
\emph{Journal of Computational and Graphical Statistics}, \bold{26} (4), 883--893.
}
\seealso{
\code{\link{depthf.fd1}}
}
\author{
Stanislav Nagy, \email{nagy at karlin.mff.cuni.cz}
}
\keyword{depth}
\keyword{functional}
\keyword{outlier}
\keyword{plot}
\keyword{shape}
