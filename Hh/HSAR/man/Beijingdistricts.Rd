\name{Beijingdistricts}
\alias{Beijingdistricts}
\docType{data}
\title{Boundaries of districts (SpatialPolygonsDataFrame) in Beijing}
\description{
  The geographic boundaries of districts (SpatialPolygonsDataFrame) in Beijing. This is only a subset of districts in Beijing where our land parcel samples fall into.

}
\usage{data(Beijingdistricts)}

\seealso{ \code{\link{landprice}}, \code{\link{landSPDF}}}

\examples{
data(Beijingdistricts)

library(spdep)
plot(Beijingdistricts,border="light grey")
# extract the area of each district
library(rgeos)
library(classInt)
library(RColorBrewer)

Beijingdistricts$geo.area <- gArea(Beijingdistricts,byid=TRUE) / 1000000
x <- Beijingdistricts$geo.area
breaks <- classIntervals(x,4,"fisher")$brks
groups <- cut(x,breaks,include.lowest=TRUE,labels=FALSE)
palette <- brewer.pal(4, "Blues")
plot(Beijingdistricts,col=palette[groups],border="grey")

# extract the district level spatial weights matrix
nb.list <- spdep::poly2nb(Beijingdistricts,queen=FALSE)
mat.list <- spdep::nb2mat(nb.list,style="W")
M <- as(mat.list,"dgCMatrix")
}
\keyword{hierarchical spatial data}
