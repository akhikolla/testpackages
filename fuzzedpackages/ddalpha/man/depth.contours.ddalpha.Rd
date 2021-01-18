\name{depth.contours.ddalpha}
\alias{depth.contours.ddalpha}
\title{
Depth Contours
}
\description{
Builds the data depth contours for multiclass 2-dimensional data using the trained classifier.
Also accessible from \code{\link{plot.ddalpha}}.
}
\usage{
depth.contours.ddalpha(ddalpha,
              main = "", xlab="", ylab = "", 
              drawplot = T, frequency=100, levels = 10, drawsep = T, ...)
}

\arguments{
  \item{ddalpha}{
DD\eqn{\alpha}-classifier (obtained by \code{\link{ddalpha.train}}).
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
  \item{drawsep}{
draws the separation on the DD-plot (currently for 2 classes and not for knn)
}
  \item{\dots}{
additional parameters passed to the depth functions and to \code{\link{plot}}
}
}

\seealso{
 \code{\link{depth.}}, 
 \code{\link{depth.contours}}, 
 \code{\link{depth.graph}}.
}
\examples{

\dontrun{

par(mfrow = c(2,2))
data(hemophilia)

ddalpha = ddalpha.train(hemophilia, depth = "none")
depth.contours.ddalpha(ddalpha, main = "data")

for (depth in c("zonoid", "Mahalanobis", "projection", "spatial")){
  ddalpha = ddalpha.train(hemophilia, depth = depth)
  depth.contours.ddalpha(ddalpha, main = depth)
}

for (depth in c("halfspace", "simplicial", "simplicialVolume")){
  ddalpha = ddalpha.train(hemophilia, depth = depth, exact = T)
  depth.contours.ddalpha(ddalpha, main = depth)
}


}
}

\keyword{ visualization }
