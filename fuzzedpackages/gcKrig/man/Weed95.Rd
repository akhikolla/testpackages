\name{Weed95}
\alias{Weed95}
\docType{data}
\title{
Counts of Weed Plants on a Field
}
\description{
The weed species Viola Arvensis was counted within circular frames each of area
0.25 square meter except for 10 missing sites in the first row,
from a 20 by 14 rectangular grid, so the total number of locations is 270.
Also, the percentages of organic matter in a soil sample are collected.
The data was studied by Christensen and Waagepetersen (2002)
to investigate whether weed occurrence could be predicted
from observations of soil texture and soil chemical properties.
}
\usage{data(Weed95)}
\format{
  A data frame with 270 observations and 6 variables.
  \describe{
    \item{\code{xloc}}{Cartesian x-coordinate of the locations (in meter).}
    \item{\code{yloc}}{Cartesian y-coordinate of the locations (in meter).}
    \item{\code{weedcount}}{Number of weed collected at the given site.}
    \item{\code{scaled Y coord}}{The scaled Y coordinate with range -1 to 1 as a
    covariate in regression.}
    \item{\code{organic}}{Another chemical component indicating the organic matter of the soil.}
    \item{\code{dummy}}{A dummy variable taking values 0 or 1. If 0 it is
     treated as observed location and 1 treated as predicted location in Christensen and Waagepetersen (2002). }
  }
}

\references{
Christensen,O. and Waagepetersen,R. (2002) Bayesian Prediction of Spatial Count Data Using
Generalized Linear Mixed Models. \emph{Biometrics, 58:280-286 }
}
\examples{
data(Weed95)
str(Weed95)
}
\keyword{Datasets}
