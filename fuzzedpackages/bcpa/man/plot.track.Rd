\name{plot.track}
\alias{plot.track}
\title{Plot Track}
\usage{
\method{plot}{track}(x, pch = 19, col = rgb(0, 0, 0, 0.2), ...)
}
\arguments{
  \item{x}{object of class "track" to plot}

  \item{pch}{default point style (filled circle)}

  \item{col}{default color (transluscent grey)}

  \item{...}{other arguments to pass to \code{plot}}
}
\description{
Default method for plotting tracks
}
\examples{
data(Simp)
is(Simp)
plot(Simp)
}

