\name{depthf.simplicialBand}
\alias{depthf.simplicialBand}

\title{
Calculate Simplicial Band Depth
}
\description{
Calculate the simplicial band depth defined by Lopez-Pintado, Sun, Lin, Genton (2014).
}
\usage{
depthf.simplicialBand(objectsf, dataf, modified = TRUE, J = NULL, 
                      range = NULL, d = 101)
}
\arguments{
  \item{objectsf}{
Functoins for which the depth should be computed; a list containing lists (functions) of two vectors of equal length, named \code{args} and \code{vals}: arguments sorted in ascending order and corresponding them values respectively.
}
  \item{dataf}{
Data sample of functoins w.r.t. which the depth should be computed; structure as for \code{objectsf}.
}
  \item{modified}{
Whether modified simplicial band depth should be computed; logical, \code{TRUE} by default.
}
  \item{J}{
How many functions to consider in each tuple of the U-statistics; integer, \code{d+1} by default.
}
  \item{range}{
The common range of the domain where the functions of objectsf and dataf are observed. Vector of length \code{2} with the left and the right end of the interval. Must contain all arguments given in objectsf and dataf.
}
  \item{d}{
Grid size to which all the functional data are transformed. For depth computation, all functional observations are first transformed into vectors of their functional values of length d corresponding to equi-spaced points in the domain given by the interval range. Functional values in these points are reconstructed using linear interpolation, and extrapolation.
}
}
\value{
A vector of depths of each of \code{objectsf} w.r.t. \code{dataf}.
}
\references{
Lopez-Pintado, Sun, Lin, Genton (2014). Simplicial band depth for multivariate data. \emph{Advances in Data Analysis and Classification} \bold{8}(3) 321--338.
}
\seealso{
\code{\link{depthf.BD}}, \code{\link{depthf.ABD}}.
}
\examples{
datafA = dataf.population()$dataf[1:20]
datafB = dataf.population()$dataf[21:50]

dataf2A = derivatives.est(datafA,deriv=c(0,1))
dataf2B = derivatives.est(datafB,deriv=c(0,1))
depthf.simplicialBand(dataf2A,dataf2B)
}
\keyword{ robust }
\keyword{ functional }
\keyword{ nonparametric }
