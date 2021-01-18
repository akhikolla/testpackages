\name{AtlanticFish}
\alias{AtlanticFish}
\docType{data}
\title{
Dataset of Mid-Atlantic Highlands Fish}
\description{
This dataset was studied by Johnson and Hoeting (2011) for analyzing pollution tolerance
in Mid-Atlantic Highlands Fish.
Pollution intolerant fish were sampled, and several stream characteristics
were measured to assess water quality at 119 sites in the Mid-Atlantic region of the United States.
All covariates of the data had been standardized to have mean 0 and variance 1.
}
\usage{data(AtlanticFish)}
\format{
  A data frame with 119 observations and 12 variables.
  \describe{
    \item{\code{LON}}{Longitude of the location.}

    \item{\code{LAT}}{Latitude of the location.}

    \item{\code{ABUND}}{Fish abundance at given locations, the discrete response.}

    \item{\code{ORDER}}{Strahler stream order, a natural covariate measuring stream size.}

    \item{\code{DISTOT}}{Watershed classified as disturbed by human activity,
                         a variable reflecting stream quality.}

    \item{\code{HAB}}{An index of fish habitat quality at the stream site, a variable reflecting stream quality.}

    \item{\code{WSA}}{Watershed area, a natural covariate.}

    \item{\code{ELEV}}{Elevation.}

    \item{\code{RD}}{Road density in the watershed, a variable reflecting stream quality.}

    \item{\code{DO}}{Concentration of dissolved oxygen in the stream at the sampling site,
                     a stream quality variable.}

    \item{\code{XFC}}{Percent of areal fish cover at the sampling site, a stream quality variable.}

    \item{\code{PCT}}{Percent of sand in streambed substrate, a stream quality variable.}
  }
}

\references{
Johnson, D. and Hoeting, J. (2011) Bayesian Multimodel Inference for Geostatistical Regression Models,
\emph{PLOS ONE, 6:e25677}.
\url{http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0025677}.
}
\examples{
data(AtlanticFish)
str(AtlanticFish)
}
\keyword{Datasets}
