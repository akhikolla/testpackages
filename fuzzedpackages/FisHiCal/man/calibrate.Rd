\name{calibrate}
\alias{calibrate}
\title{
Hi-C calibration
}
\description{
A function for applying a calibration function to a given Hi-C contact matrix.
}
\usage{
calibrate(hic, calib)
}
\arguments{
  \item{hic}{
A numeric matrix/vector of pairwise Hi-C contact frequencies. 
}
  \item{calib}{
  A list with f, the calibration function, and params a list of parameters for f.
  This object can be generated with the function \code{\link{prepareCalib}}. 
}
}

\value{
A calibrated version of the input matrix (or vector), i.e. the result of applying the calibration function (f) on the values of 'hic'. Zero values represent discarded/missing information.
}
\references{
Y. Shavit, F.K. Hamey, P. Lio', FisHiCal: an R package for iterative FISH-based calibration of Hi-C data, 2014 (submitted).
}
\author{
Yoli Shavit
}

\seealso{
\code{\link{prepareCalib}}\cr
\code{\link{conf}}\cr
\code{\link{hic}}\cr
\code{\link{match}}
}

\examples{
 data(match) 
 data(hic)
 npoints = 10 # number of points to fit
 
 # prepareCalib computes threshold to fit 
 res = prepareCalib(match, npoints, useMax = FALSE)
 calib = res$calib
 calib$f # calibration function
 calib$params # the parametrs for f
 # note that calib could be refined by the user
 
 # now calibrate the hic matrix
 calibHiC = calibrate(hic, calib)
 plot(match$distances, calibHiC[upper.tri(calibHiC)],
 xlab = "distances", ylab = "calibrated distances")
}
