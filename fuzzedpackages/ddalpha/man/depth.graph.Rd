\name{depth.graph}
\alias{depth.graph}
\title{
Depth Graph
}
\description{
Builds the data depth graphs for 2-dimensional data. The graph is built using \code{\link{persp}}.
}
\usage{
depth.graph(data, 
  depth_f = c("halfspace", "Mahalanobis", "projection", "simplicial", 
              "simplicialVolume", "spatial", "zonoid", "none"), 
  apoint = NULL, 
  main = depth_f,  
  xlim = c(min(data[, 1]), max(data[, 1])), 
  ylim = c(min(data[, 2]), max(data[, 2])), 
  zlim = c(0, max(z)), 
  xnum = 250, 
  ynum = 250, 
  theta=15, phi=60,
  bold = F,
  ...)
}

\arguments{
  \item{data}{
2-dimensional numeric data frame or matrix
}
  \item{depth_f}{
the name of the depth function. The list of the supported depths and described in the topic \code{\link{depth.}}.
}
  \item{apoint}{
a 2-dimensional point which is shown in black color.
}
  \item{main}{
an overall title for the plot: see \code{\link{title}}
}
  \item{xlim, ylim, zlim}{
numeric vectors of length 2, giving the x, y and z coordinates ranges: see \code{\link{plot.window}}
}
  \item{xnum, ynum}{
number of points on each direction, x and y. Impacts the smoothness of the surface. 
}
  \item{theta, phi}{
rotation angles
}
  \item{bold}{
draws bold points
}
  \item{\dots}{
additional parameters passed to \code{\link{persp}}
}
}

\seealso{
\code{\link{depth.}}

\code{\link{persp}}
}
\examples{

\dontrun{

par(mfrow = c(2,3), mar = c(0,0,0,0), mai = c(0,0,0.2,0))
data(hemophilia)
depth.graph(hemophilia, "none", xnum = 100, ynum = 100)
depth.graph(hemophilia, "Mahalanobis", xnum = 100, ynum = 100)
depth.graph(hemophilia, "halfspace", xnum = 100, ynum = 100)
depth.graph(hemophilia, "projection", xnum = 100, ynum = 100)
depth.graph(hemophilia, "zonoid", xnum = 100, ynum = 100)
depth.graph(hemophilia, "spatial", xnum = 100, ynum = 100)

}
}

\keyword{ visualization }
