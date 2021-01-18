\name{plot.ddalphaf}
\alias{plot.ddalphaf}
\title{
Plots for the "ddalphaf" Class
}
\description{

\code{\link{plot.functional}} -- plots the functional data used by classifier

\code{\link{depth.contours.ddalpha}} -- builds the data depth contours for multiclass 2-dimensional data using the trained classifier.
\code{\link{draw.ddplot}} -- draws the \emph{DD}-plot of the existing DD\eqn{\alpha}-classifier.
}
\usage{
\method{plot}{ddalphaf}(x, type = c("functional.data", "ddplot", "depth.contours"), ...)
}

\arguments{
  \item{x}{
functional DD\eqn{\alpha}-classifier (obtained by \code{\link{ddalphaf.train}}).
}
  \item{type}{
type of the plot for \code{\link{plot.functional}}, \code{\link{draw.ddplot}} or \code{\link{depth.contours.ddalpha}}
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

dataf = dataf.growth()
ddalphaf = ddalphaf.train (dataf$dataf, dataf$labels, 
                            classifier.type = "ddalpha", maxNumIntervals = 2)

# plot the functional data
plot(ddalphaf)

# plot depth contours and separation in the transformed space 
# (possible only if maxNumIntervals = 2)
plot(ddalphaf, type = "depth.contours")

# plot the DD-plot
plot(ddalphaf, type = "ddplot")

}
}
\keyword{ functional }
\keyword{ visualization }
