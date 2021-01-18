\name{sbmSDP-package}
\alias{sbmSDP-package}
\alias{sbmSDP}
\docType{package}
\title{
Semidefinite Programming for Fitting Block Models of Equal Block Sizes
}
\description{
An ADMM implementation of SDP-1, a semidefinite programming relaxation of 
the maximum likelihood estimator for fitting a block model. SDP-1 has a tendency 
to produce equal-sized blocks and is ideal for producing a form of network histogram
approximating a nonparametric graphon model. Alternatively, it can be used for community detection. 
(This is experimental code, proceed with caution.)
}
\details{
\tabular{ll}{
Package: \tab sbmSDP\cr
Type: \tab Package\cr
Version: \tab 0.2\cr
Date: \tab 2015-06-18\cr
License: \tab GPL-3\cr
}
An ADMM implementation of SDP-1 algorithm for fitting stochastic block models (SBMs). The main function is sdp1_admm.
}
\author{
Arash A. Amini

Maintainer: Arash A. Amini <amini.aa@gmail.com>
}
\references{
  On Semidefinite relaxations of the block model by A.A. Amini and E. Levina.
}
\keyword{ 
community detection
stochastic block model 
}
