\name{OilWell}
\alias{OilWell}
\docType{data}
\title{
Location of Successful and Dry Wells
}
\description{
A dataset recording locations of successful and unsuccessful drilling oil wells
in the northwest shelf of Delaware basin in New Mexico, a region that is
densely drilled but has some sparsely drilled areas.
The original dataset was transformed to a central area of about 65 square kilometers,
see Hohn (1999), Chapter 6.
}
\usage{data(OilWell)}
\format{
  A data frame with 333 observations and 3 variables.
  \describe{
    \item{\code{Easting}}{Cartesian x-coordinate of the locations.}
    \item{\code{Northing}}{Cartesian y-coordinate of the locations.}
    \item{\code{Success}}{A binary variable indicating the success of the drill at given locations.
                        1 for successful drilling oil wells and 0 for unsuccess. }
  }
}

\references{
Hohn, M. (1999) Geostatistics and petroleum geology, Second Edition.
\emph{Springer Science & Business Media, APA}
}
\examples{
data(OilWell)
str(OilWell)
}
\keyword{Datasets}
