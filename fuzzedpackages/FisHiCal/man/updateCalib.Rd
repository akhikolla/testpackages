\name{updateCalib}
\alias{updateCalib}
\title{
Update a parameter value in the calibration model
}
\description{
A function to update a given parameter in the calibration model. 
}
\usage{
updateCalib(calib, paramVal, paramIndex)
}
\arguments{
\item{calib}{
  A list with f, the calibration function, and params a list of parameters for f.
  This object can be generated with the function \code{\link{prepareCalib}}. } 
\item{paramVal}{
The new value of the parameter to be updated.
}
\item{paramIndex}{
The index of the parameter to be updated in the parameter list of calib.
}
}

\details{
Users may want to refine the parameters of the calibration model. For example, in the model generated with \code{\link{prepareCalib}}\cr the user may wish to increase the threshold value in order to discard less information in the calibration step. 
}                                                                   

\value{
An updated calibration model (an updated calib object).
}

\references{
Y. Shavit, F.K. Hamey, P. Lio', FisHiCal: an R package for iterative FISH-based calibration of Hi-C data, 2014 (submitted).
}
\author{
Yoli Shavit
}
\seealso{
\code{\link{prepareCalib}}\cr
\code{\link{calibrate}}\cr
\code{\link{match}}\cr
}

\examples{ 
 data(match) 
 npoints = 10 # number of points to fit
 
 # prepareCalib computes threshold according to the fit 
 # useMax is set to FALSE
 res = prepareCalib(match, npoints, useMax = FALSE)
 calib = res$calib
 calib
 calib = updateCalib(calib, 3.9, 3)
 calib
             
}


