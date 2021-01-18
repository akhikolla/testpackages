\docType{package}

\name{psgp}

\alias{psgp-package}

\title{Projected Sequential Gaussian Processes for spatial interpolation}

\author{Ben Ingram, Remi Barillec}

\description{ The \code{psgp} package provides a spatial interpolation method based on Projected Sequential Gaussian Processes (PSGP) for the \code{intamap} package. The PSGP (Csato, 2002) is an approximation to the standard Gaussian process (Rasmussen and Williams, 1996) whereby the observations are projected onto a subset of optimal "active" observations, thus reducing possible redundancy in the data and allowing for faster, memory efficient, interpolation. The projection is done in a sequential manner, that is one observation is projected onto the active subset at a time. This allows for larger datasets to be processed, and overcomes the limitations of standard Gaussian processes related to storing the full covariance matrix (which can be unfeasible for really large datasets).

This implementation of PSGP for spatial interpolation uses a mixture of covariance kernels, namely an Exponential kernel and a Matern kernel with roughness parameter 5/2. The covariance function also includes a nugget term (white noise) and bias term.
}

\section{Parameter estimation}{ The method \code{\link{estimateParameters}} looks for 
  a set of covariance function parameters (kernels, white noise and bias) which maximise
  the likelihood of the observation.}

\section{Prediction/Interpolation}{ The method \code{\link{spatialPredict}} computes predictions (including variance) at a set of unobserved locations.}

\section{Measurement errors}{ It is possible to include measurement errors if these are available. These will be taken into account when estimating parameters and making predictions. See \code{\link{learnParameters}} for more information. }

\section{System requirements}{ the PSGP package is written in C++ and uses the Armadillo library for all linear algebra routines.}
\references{
  
  Csato, L., Gaussian Processes - Iterative Sparse Approximations, Ph.D Thesis, NCRG, Aston University, UK, 2002.
  
  Rasmussen, C. E. and Williams, C. K., Gaussian Processes for Machine Learning,  The MIT Press, Cambridge, Massachusetts, 1996.
  
  \url{http://arma.sourceforge.net}
}
