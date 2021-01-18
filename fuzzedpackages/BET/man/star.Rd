\docType{data}
\name{star}
\alias{star}
\title{Coordinates of Brightest Stars in the Night Sky}
\format{
An object of class \code{data.frame} with 256 rows and 2 columns.
}
\usage{
data(star)
}
\description{
This data set collects the galactic coordinates of the 256 brightest stars in the night sky (Perryman et al. 1997). We consider the longitude (\code{x}) and sine latitude (\code{y}) here.
}
\examples{
data(star)
BETs(cbind(star$x.raw, star$y.raw), asymptotic = FALSE, plot = TRUE)
}
\keyword{datasets}

