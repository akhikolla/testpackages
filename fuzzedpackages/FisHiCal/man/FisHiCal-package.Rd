\name{FisHiCal-package}
\alias{FisHiCal-package}
\alias{FisHiCal}
\docType{package}
\title{
Iterative FISH-based Calibration of Hi-C Data
}
\description{
FisHiCal is an R package for combining Hi-C and FISH data, which offers a modular and easy-to-use tool for chromosomal spatial analysis. With FisHiCal researchers can prepare and apply FISH-based Hi-C calibration, which converts contact frequencies into distances and computes a range threshold, and exploit 3D inference methods to iteratively refine it. These methods include algorithms for reconstructing chromosome structure (from calibrated distances) and for detecting spatial inconsistencies. 
}
\details{
\tabular{ll}{
Package: \tab FisHiCal\cr
Type: \tab Package\cr
Version: \tab 1.1\cr
Date: \tab 2014-06-20\cr
License: \tab GPL \cr
}
FisHiCal implements the methods described in Shavit, Hamey and Lio' (2014) (submitted). Users can first prepare and apply their calibration with \code{\link{prepareData}}, \code{\link{prepareCalib}} and \code{\link{calibrate}} and use 3D inference functions (\code{\link{lsmacof}}, \code{\link{searchInc}}, \code{\link{plotInc}} and \code{\link{summaryInc}}) to spatially explore their data and further refine their calibration.
}
\author{
Yoli Shavit, Fiona Kathryn Hamey and Pietro Lio'.\cr
Maintainer: Yoli Shavit <ys388@cam.ac.uk>
}
\references{
Y. Shavit, F.K. Hamey, P. Lio', FisHiCal: an R package for iterative FISH-based calibration of Hi-C data, 2014 (submitted).
}
\keyword{ package }

