\name{LansingTrees}
\alias{LansingTrees}
\docType{data}
\title{
Locations and Botanical Classification of Trees in Lansing Woods
}
\description{
The data is aggregated from the dataset \code{lansing} in library \code{spatstat}, which
came from an investigation of a 924 ft x 924 ft (19.6 acres) plot in Lansing Woods, Clinton
County, Michigan USA by D.J. Gerrard.
The original point process data described the locations of 2,251 trees and their
botanical classification (into maples, hickories, black oaks, red oaks, white oaks and miscellaneous trees). The original plot size has been rescaled to the unit square and the number of different types of trees
has been counted within squares of length 1/16.}

\usage{data(LansingTrees)}
\format{
  A data frame with 256 observations and 8 variables.
  \describe{
    \item{\code{Easting}}{Cartesian x-coordinate of the locations.}

    \item{\code{Northing}}{Cartesian y-coordinate of the locations.}

    \item{\code{maple}}{Number of maples in the area.}
    \item{\code{hickory}}{Number of hickories in the area.}
    \item{\code{blackoak}}{Number of black oaks in the area.}
    \item{\code{redoak}}{Number of red oaks in the area.}
    \item{\code{whiteoak}}{Number of white oaks in the area.}
    \item{\code{misc}}{Number of miscellaneous trees in the area.}
  }
}


\references{
Kazianka, H. (2013) Approximate copula-based estimation and prediction of discrete spatial data.
\emph{Stoch Environ Res Risk Assess, 27:2015-2026}
}
\examples{
data(LansingTrees)
str(LansingTrees)
}
\keyword{Datasets}
