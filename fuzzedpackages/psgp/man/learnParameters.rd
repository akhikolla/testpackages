\name{learnParameters}
\alias{learnParameters}
\title{Projected Sequential Gaussian Process}
\description{
  \code{learnParameters} performs maximum likelihood parameter estimation in the PSGP framework.
}
\usage{
  learnParameters(object)
}
\arguments{
  \item{object}{ a list object of intamap type. Most arguments necessary for
    interpolation are passed through this object. 
    See \code{\link[intamap:intamap-package]{intamap-package}} for further description of the necessary 
    content of this variable.
  }
}

\section{Warning}{It is advised to use the intamap wrapper \code{\link{estimateParameters}} rather than calling this method directly.}

\details{
  The Projected Spatial Gaussian Process (PSGP) framework provides  
  an approximation to the full Gaussian process in which the observations 
  are projected sequentially onto an optimal subset of 'active' observations. Spatial 
  interpolation is done using a mixture of covariance kernels (exponential and Matern 
  5/2).

  The function \code{learnParameters} is an internal function for estimating the
  parameters of the covariance function given the data, using a maximum likelihood
  approach. A valid intamap \code{object} must be passed in.
  
  PSGP is able to also take the measurement characteristics (i.e. errors) into
  account using possibly many error models. For each error model, assumed Gaussian, the 
  error variance can be specified. The vector 
  \code{object$observations$oevar} contains all variances for the error models (one 
  value per error model). 
  % Similarly, if specified, the vector  
  %\code{object$observations$oebias} must contains the bias values for each model. 
  Which error model is used for a given observation is determined by the 
  \code{object$observations$oeid} vector of indices, which specifies the index of the 
  model to be used for each observation.
}


\references{ 
    Csato and Opper, 2002; Ingram et al., 2008
    
}
\author{Ben Ingram, Remi Barillec}
\seealso{
  \code{\link{makePrediction}},
  \code{\link{learnParameters}},
  \code{\link[intamap:estimateParameters]{estimateParameters}},
  \code{\link[intamap:createIntamapObject]{createIntamapObject}}
}
\examples{
  # see example in estimateParameters
}
\keyword{spatial}


