\name{MakeTrack}
\alias{MakeTrack}
\title{Make Track}
\usage{
MakeTrack(X, Y, Time)
}
\arguments{
  \item{X}{vector of X locations}

  \item{Y}{vector of Y locations}

  \item{Time}{vector of time (can be POSIX)}
}
\value{
a \code{track} class data frame, with three columns:
\code{X}, \code{Y} and \code{Time}.
}
\description{
Simple convenience function for creating a \code{track}
class object from X-Y-Time movement data. A \code{track}
class object can be conveniently plotted and analyzed
within \code{bcpa}.
}
\examples{
X <- cumsum(arima.sim(n=100, model=list(ar=0.8)))
Y <- cumsum(arima.sim(n=100, model=list(ar=0.8)))
Time <- 1:100
mytrack <- MakeTrack(X,Y,Time)
plot(mytrack)
}
\seealso{
plot.track
}

