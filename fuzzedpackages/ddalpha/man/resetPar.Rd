\name{resetPar}
\alias{resetPar}
\title{
Reset Graphical Parameters
}
\description{
The function returns the default graphical parameters for \code{\link{par}}.
}
\usage{
resetPar()
}
\details{
The returned parameters are used as input parameters for \code{\link{par}}.
}
\value{
The list of graphical parameters described in the 'Graphical Parameters' section of \code{\link{par}}.
}
\examples{

par(mfrow = c(1,2), mar = c(2,2,2,2))
plot(sin, -pi, 2*pi)
plot(cos, -pi, 2*pi)

par(resetPar())
plot(sin, -pi, 2*pi)
plot(cos, -pi, 2*pi)

}

\keyword{ visualization }
