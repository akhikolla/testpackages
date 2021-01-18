\name{buildMetadata}
\alias{buildMetadata}
\title{Build Metadata Table}
\description{
  \code{buildMetadata} builds a metadata table of likelihood model descriptors in the PSGP framework. This is an internal function and it should not be used directly.
}
\usage{
  buildMetadata(observations)
}
\arguments{
  \item{observations}{ an observations data frame containing a vector of observation variances (\code{obserations$oevar}) and, optionally, a list of observation biases (\code{observationS$oebias}). If the biases are omitted, a zero bias is assumed for all likelihood models.}
}
\details{
  \code{buildMetadata} builds a metadata table of likelihood model descriptors in the PSGP framework. The likelihood models are assumed Gaussian with variances specified in the vector \code{observations$oevar} (the bias is assumed to be zero). Optionally, biases can be specified in the \code{observations$oebias} vector. However, biases are not taken into account in the current version of the PSGP package (they will be in a future release). 
}

\author{Remi Barillec}
\seealso{
  \code{\link{learnParameters}},
  \code{\link{makePrediction}},
  \code{\link[intamap:estimateParameters]{estimateParameters}},
  \code{\link[intamap:spatialPredict]{spatialPredict}},
}
\examples{
  ## Load our favourite dataset
  data(meuse)
  obs <- meuse
  
  ## Number of observations
  nobs <- length(obs$y)
  
  ## Indicate which likelihood model should be used for each observation
  obs$oeid  <- seq(1:nobs)
  
  ## Use random variances for the sake of the example
  obs$oevar <- rnorm( max(obs$oeid) )
  
  ## Generate metadata table and print it out
  metadata <- buildMetadata(obs)
  print(metadata)
}
\keyword{spatial}



