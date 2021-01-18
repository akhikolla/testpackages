\name{makePrediction}
\alias{makePrediction}
\title{Spatial projected sequential GP prediction}
\description{\code{makePrediction} performs prediction/interpolation within the PSGP package.
}
\usage{ makePrediction(object, vario) }
\arguments{
  \item{object}{ a list object of intamap type. Most arguments necessary for
    interpolation are passed through this object. See \link[intamap:intamap-package]{intamap-package} for 
    further description of the necessary content of this variable.  Additional meta data 
    about the measurement process is included in this object. In particular, see 
    \link{learnParameters} for a way to specify measurement error variances.
  }
  \item{vario}{ Log-parameters of the covariance function. For compatibility with the 
    intamap package, the log-parameters of the PSGP covariance function are stored within
    a variogram array object (see \code{\link[gstat:vgm]{vgm}}), as follows: \cr
    vario[1,1]   NA \cr
    vario[1,2]   length scale (or range) of the Exponential kernel \cr
    vario[1,3]   variance (or sill) of the Exponential kernel \cr
    vario[1,4]   length scale (or range) of the Matern 5/2 kernel \cr
    vario[1,5]   variance (or sill) of the Matern 5/2 kernel \cr
    vario[1,6]   inverse bias (i.e. 1/mean(data)) \cr
    vario[1,7]   white noise variance (nugget)
  } 
} 

\section{Warning}{It is advised to use the intamap wrapper \code{\link{spatialPredict}} rather than calling this method directly.}

\details{
  The Projected Spatial Gaussian Process (PSGP) framework provides  
  an approximation to the full Gaussian process in which the observations 
  are projected sequentially onto an optimal subset of 'active' observations. Spatial 
  interpolation is done using a mixture of covariance kernels (Exponential and Matern 5/2).
  
  The function \code{makePrediction} is a function for making predictions at a set 
  of unobserved inputs (or locations).  
  
  Measurement characteristics (i.e. observation error) can be specified if needed. 
  See \link{learnParameters} for a description of how to specify measurement error
  models with given variances.
}

\references{ 
 L. Csato and M. Opper. Sparse online Gaussian processes. Neural Computation, 14(3):
641-669, 2002.

 B. Ingram, D. Cornford, and D. Evans. Fast algorithms for automatic mapping with space-
limited covariance functions. Stochastic Environmental Research and Risk Assessment, 22
(5):661-670, 2008.
}
\author{Ben Ingram, Remi Barillec}
\seealso{
  \code{\link{learnParameters}}
  \code{\link[intamap]{spatialPredict}}, 
  \code{\link[intamap]{createIntamapObject}}
}
\examples{
  # see example in spatialPredict
}
\keyword{spatial}
