\name{landSPDF}
\alias{landSPDF}
\docType{data}
\title{Beijing land price data}
\description{
  The spatial locations (SpatialPointsDataFrame) of the Beijing land price data. It is used to extract the spatial weights matrix among land parcels.   

}
\usage{data(landSPDF)}

\seealso{ \code{\link{landprice}}, \code{\link{Beijingdistricts}}}

\examples{
data(landSPDF)

library(spdep)
plot(landSPDF,col="red",pch=16,cex=0.7)

data(landprice)

library(classInt)
library(RColorBrewer)

# link the variables in the landprice to the spatial data
index.match <- match(landSPDF$obs,landprice$obs)
landSPDF@data <- data.frame(landSPDF@data,landprice[index.match,])

par(mar=c(0,0,0,0))
x <- landSPDF$lnprice
breaks <- classIntervals(x,4,"fisher")$brks
groups <- cut(x,breaks,include.lowest=TRUE,labels=FALSE)
palette <- brewer.pal(4, "Blues")
plot(landSPDF,pch=19,col=palette[groups],cex=0.8)

# extract a spatial weights matrix based on the distances between pairs of land parcels
nb.25 <- spdep::dnearneigh(landSPDF,0,2500)
# to a weights matrix
dist.25 <- spdep::nbdists(nb.25,landSPDF)
dist.25 <- lapply(dist.25,function(x) exp(-0.5 * (x / 2500)^2))
mat.25 <- spdep::nb2mat(nb.25,glist=dist.25,style="W")
W <- as(mat.25,"dgCMatrix")
}
\keyword{hierarchical spatial data}
