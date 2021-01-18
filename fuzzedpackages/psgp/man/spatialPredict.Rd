\name{spatialPredict.psgp}
\alias{spatialPredict.psgp}
\alias{spatialPredict}
\title{Spatial prediction using a Projected Sequential Gaussian Process (PSGP)}
\description{This overloads the \code{\link[intamap:spatialPredict]{spatialPredict}} routine 
  from the intamap package for interpolation using the PSGP method.
}
\usage{ spatialPredict(object, ...) }
\arguments{
  \item{object}{ a list object of type PSGP. Most arguments necessary for 
    interpolation are passed through this object. See \code{\link[intamap:intamap-package]{intamap-package}} 
    for further description of the necessary content of this variable
  }
  \item{...}{ optional extra arguments (these are only used for debugging purposes)}
} 

\details{
  See \code{\link{psgp-package}} and \code{\link{makePrediction}} for 
  further detail.
}

\author{Ben Ingram, Remi Barillec}
\seealso{
  \code{\link{psgp-package}}, 
  \code{\link[intamap:estimateParameters]{estimateParameters}}, 
  \code{\link{makePrediction}}
  \code{\link[intamap:createIntamapObject]{createIntamapObject}}
}

\examples{
data(meuse)
meuse = meuse[1:100,]
coordinates(meuse) = ~x+y
meuse$value = log(meuse$zinc)
data(meuse.grid)
gridded(meuse.grid) = ~x+y
proj4string(meuse) = CRS("+init=epsg:28992")
proj4string(meuse.grid) = CRS("+init=epsg:28992")

# Specify a different observation error model for each observation 
nobs = length(meuse$value)      # Number of observations
meuse$oeid = seq(1:nobs)        # One error model per observation

# Indicate the variance for each of these error models
meuse$oevar <- abs( rnorm( max(meuse$oeid) ) )

# Set up intamap object
obj = createIntamapObject(
  observations = meuse,
  predictionLocations = meuse.grid,
  targetCRS = "+init=epsg:3035",
  class = "psgp"
)

# Estimate parameters and predict at new locations (interpolation)
obj = conformProjections(obj)
obj = estimateParameters(obj)
obj = spatialPredict(obj)

# Plot results
plotIntamap(obj)
}
\keyword{spatial}
