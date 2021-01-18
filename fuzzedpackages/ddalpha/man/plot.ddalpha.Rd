\name{plot.ddalpha}
\alias{plot.ddalpha}
\title{
Plots for the "ddalpha" Class
}
\description{
\code{\link{depth.contours.ddalpha}} -- builds the data depth contours for multiclass 2-dimensional data using the trained classifier.
\code{\link{draw.ddplot}} -- draws the \emph{DD}-plot of the existing DD\eqn{\alpha}-classifier.
}
\usage{
\method{plot}{ddalpha}(x, type = c("ddplot", "depth.contours"), ...)
}

\arguments{
  \item{x}{
DD\eqn{\alpha}-classifier (obtained by \code{\link{ddalpha.train}}).
}
  \item{type}{
type of the plot for \code{\link{draw.ddplot}} or \code{\link{depth.contours.ddalpha}}
}
  \item{\dots}{
additional parameters passed to the depth functions and to \code{\link{plot}}
}
}

\seealso{
\code{\link{depth.}}

\code{\link{depth.contours}}

\code{\link{depth.graph}}
}
\examples{

\dontrun{


par(mfrow = c(2,2))
data(hemophilia)

ddalpha = ddalpha.train(hemophilia, depth = "none")
plot(ddalpha, type = "depth.contours", main = "data")
plot(ddalpha, type = "ddplot", main = "data", drawsep = F)

for (depth in c("zonoid", "Mahalanobis", "projection", "spatial")){
  ddalpha = ddalpha.train(hemophilia, depth = depth)
  plot(ddalpha, type = "depth.contours", main = depth, drawsep = T)
  plot(ddalpha, type = "ddplot", main = depth)
}

}
}

\keyword{ visualization }
