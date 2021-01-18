\name{estimateParameters.psgp}
\alias{estimateParameters.psgp}
\alias{estimateParameters}
\title{Parameter estimation using a Projected Sequential Gaussian Process (PSGP)}
\description{
This overloads the \code{\link[intamap:estimateParameters]{estimateParameters}} routine 
  from the intamap package for interpolation using the PSGP method.
}
\usage{
  estimateParameters(object, ...)
}
\arguments{
  \item{object}{ a list object of Intamap type. Most arguments necessary for
    interpolation are passed through this object. 
    See \link[intamap:intamap-package]{intamap-package} for further description of the necessary 
    content of this variable.
  }
  \item{...}{ other parameters for the generic method, not used for this method}
}

\details{
  See \code{\link{psgp-package}} and \code{\link{learnParameters}} for 
  further details.
}

\author{Remi Barillec, Ben Ingram}
\seealso{
  \code{\link{learnParameters}},
  \code{\link[intamap:estimateParameters]{estimateParameters}},
  \code{\link{makePrediction}},
  \code{\link[intamap:createIntamapObject]{createIntamapObject}}
}
\examples{
# load our favourite dataset
data(meuse)
coordinates(meuse) = ~x+y
meuse$value = log(meuse$zinc)
data(meuse.grid)
gridded(meuse.grid) = ~x+y
proj4string(meuse) = CRS("+init=epsg:28992")
proj4string(meuse.grid) = CRS("+init=epsg:28992")

# the following two steps are only needed if one wishes to
# include observation errors

# indicate which likelihood model should be used for each observation
# in this case we use a different model for each observation
nobs = length(meuse$value)          # Number of observations
meuse$oeid  <- seq(1:nobs)
  
# the variances for the error models are random in this example
# in real examples they will come from actual measurements 
# characteristics
meuse$oevar <- abs( rnorm( max(meuse$oeid) ) )

# set up intamap object:
obj = createIntamapObject(
  observations = meuse,
  predictionLocations = meuse.grid,
  targetCRS = "+init=epsg:3035",
  class = "psgp"    # Use PSGP for parameter estimation/interpolation
)

# do interpolation step:
obj = conformProjections(obj)
obj = estimateParameters(obj)
}
\keyword{spatial}


