\name{landscapeTEST}
\alias{landscapeTEST}
\alias{landscapeTEST1}
\alias{landscapeTEST2}
\alias{landscapeTEST3}
\alias{landscapeTEST4}
\alias{landscapeTEST5}
\docType{data}
\title{
  Landscapes
}
\description{
  Five simulated landscapes, composed of 155, 154, 152, 153 and 156 fields, respectively.
}
\usage{
 landscapeTEST1
 landscapeTEST2
 landscapeTEST3
 landscapeTEST4
 landscapeTEST5}
\format{
  Landscapes have been generated using a T-tesselation algorithm. 
  The format is a formal class 'SpatialPolygons' [package "sp"].
}
\details{
  The landscape structure is simulated using a T-tessellation algorithm (Kiêu et al. 2013) in order to control specific features such as number, area and shape of the fields.
}
\references{
  Kiêu K, Adamczyk-Chauvat K, Monod H, Stoica RS. A completely random T-tessellation model and Gibbsian extensions. Spat. Stat. 2013;6:118-38. doi: http://dx.doi.org/10.1016/j.spasta.2013.09.003.
}
\examples{
  library(sp)
  library(landsepi)
  landscapeTEST1
  plot(landscapeTEST1)
}
\keyword{datasets}
