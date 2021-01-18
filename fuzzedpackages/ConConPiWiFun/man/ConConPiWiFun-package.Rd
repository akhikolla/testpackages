\name{ConConPiWiFun-package}
\alias{ConConPiWiFun-package}
\alias{ConConPiWiFun}
\docType{package}
\title{
This package contains an implementation of continuous convex piecewise (linear) functions (quadratic coming soon)
}
\description{
Continuous convex piecewise linear (ccpl) resp. quadratic (ccpq) functions can be implemented with sorted breakpoints and slopes. This includes functions that are ccpl (resp. ccpq) on a convex set (i.e. an interval or a point) and infinite out of the domain. These functions can be very usefull for a large class of optimisation problems. Efficient manipulation (such as log(N) insertion) of such data structure is obtained with map standard template library of C++ (that hides balanced trees). This package is a wrapper on such a class based on Rcpp modules. 
}
\details{
\tabular{ll}{
Package: \tab ConConPiWiFun\cr
Type: \tab Package\cr
Version: \tab 0.3.0\cr
Date: \tab 2013-02-08\cr
License: \tab GPL\cr
}
}
\author{
Robin Girard

Maintainer: <robin.girard@mines-paristech.fr>
}
\references{
Related Papers are  
}
\keyword{Optimisation, Dynamic programming}

\examples{
library(ConConPiWiFun)
####  See 
#? cplfunction for continuous convex piecewise functions
#? cplfunctionvec for (optimized) list of continnuous convex piecewise functions



}
