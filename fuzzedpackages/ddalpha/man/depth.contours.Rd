\name{depth.contours}
\alias{depth.contours}
\title{
Depth Contours
}
\description{
Builds the data depth contours for 2-dimensional data.
}
\usage{
depth.contours(data, depth, 
              main = "", xlab="", ylab = "", 
              drawplot = T, frequency=100, levels = 10,
              col = "red",
              ...)
}

\arguments{
  \item{data}{
2-dimensional numeric data frame or matrix
}
  \item{depth}{
the name of the depth function. The list of the supported depths and described in the topic \code{\link{depth.}}.
}
  \item{main}{
an overall title for the plot: see \code{\link{title}}
}
  \item{xlab, ylab}{
labels of the axes
}
  \item{drawplot}{
if set to false, the contours are built on the existing plot. 
}
  \item{frequency}{
number of points on each direction, x and y. Impacts the smoothness of the contours.
}
  \item{levels}{
numeric vector of levels at which to draw contour lines. 
If the vector contains only ONE element, the levels are generated automatically as \code{seq(0, max(depth), length.out = levels)}.
}
  \item{col}{
color, used to draw points and contours
}
  \item{\dots}{
additional parameters passed to the depth functions and to \code{\link{plot}}
}
}

\seealso{
 \code{\link{depth.}}, 
 \code{\link{depth.contours.ddalpha}}, 
 \code{\link{depth.graph}}.
}
\examples{

\dontrun{

par(mfrow = c(2,2))
data(hemophilia)

depth.contours(hemophilia[,1:2], depth = "none", main = "data")

for (depth in c("zonoid", "Mahalanobis", "projection", "spatial")){
  depth.contours(hemophilia[,1:2], depth = depth, main = depth)
}

for (depth in c("halfspace", "simplicial", "simplicialVolume")){
  depth.contours(hemophilia[,1:2], depth = depth, main = depth, exact = T)
}


}
}

\keyword{ visualization }
