\name{SpatialBSS-package}
\alias{SpatialBSS-package}
\docType{package}
\title{
Blind Source Separation for Multivariate Spatial Data
}

\description{
Blind source separation for multivariate spatial data based on simultaneous/joint diagonalization of local covariance matrices. This package is an implementation of the methods described in Nordhausen, Oja, Filzmoser, Reimann (2015) <https://doi.org/10.1007/s11004-014-9559-5> and Bachoc, Genton, Nordhausen, Ruiz-Gazen and Virta (2020) <doi:10.1093/biomet/asz079>.
}

\details{
\tabular{ll}{
Package: \tab SpatialBSS\cr
Type: \tab Package\cr
Version: \tab 0.9-0 \cr
Date: \tab 2020-11-04\cr
License: \tab GPL (>= 2)\cr
}

This package provides functions to solve the Blind Source Separation problem for multivariate spatial data. These methods are designed to work with random fields that are observed on irregular locations. Moreover, the random field is assumed to show weak second order stationarity. The main function of this package is:
\itemize{
\item{\code{\link{sbss}}} {This function derives a set of local scatter matrices that are based on spatial kernel functions, where the spatial kernel functions can be chosen. Then this set of local covariance matrices as well as the sample covariance matrix are simultaneously/jointly diagonalized. Local covariance matrices as well as local difference matrices are implemented.}
}

Joint diagonalization is computed with the \code{\link[JADE]{frjd}} (fast real joint diagonalization) algorithm from the package \code{\link[JADE:JADE-package]{JADE}}.

The random field can be either a pair of numeric matrices giving the coordinates and field values or an object of class \code{\link[sp]{SpatialPointsDataFrame}} or \code{\link[sf]{sf}}.

}

\author{
Christoph Muehlmann, Klaus Nordhausen, Joni Virta

Maintainer: 
Christoph Muehlmann \email{christoph.muehlmann@tuwien.ac.at}
}

\references{
Bachoc, F., Genton, M. G, Nordhausen, K., Ruiz-Gazen, A. and Virta, J. (2020), \emph{Spatial Blind Source Separation}, Biometrika, \url{https://doi.org/10.1093/biomet/asz079}, \cr <doi:10.1093/biomet/asz079>.

Nordhausen, K., Oja, H., Filzmoser, P., Reimann, C. (2015), \emph{Blind Source Separation for Spatial Compositional Data}, Mathematical Geosciences 47, 753-770, \url{https://doi.org/10.1007/s11004-014-9559-5}, <https://doi.org/10.1007/s11004-014-9559-5>.

}

\keyword{ package }
