\name{prepareCalib}
\alias{prepareCalib}
\title{
Prepare Hi-C calibration 
}
\description{
A function to build a calibration function, by fitting a subset of FISH distances and Hi-C frequencies with a power law model (see details). The number of distances to fit (taking distances by increasing order) or a subset of selected distances should be provided by the user. Users can also choose how to estimate the distance threshold or may explicitly provide one. 
}
\usage{
prepareCalib(data, npoints, threshold = NULL, useMax = TRUE, delta = 0.05, buffer = 1.0)
}
\arguments{
  \item{data}{
A data frame with 2 mandatory columns: distances and frequencies, standing for matching FISH distances and Hi-C frequencies, correspondingly.
This data structure could be prepared with  \code{\link{prepareData}}
}
  \item{npoints}{
An integer or an integer vector. If an integer n is given, than the shortest n distances and their matching frequencies will be used. Otherwise, the indices in the integer vector will indicate the subset of distances and frequencies to use from 'data'.   
}
  \item{threshold}{
Optional numeric, set to NULL by default. If provided, will be used as the distance range threshold of the calibration. 
}
\item{useMax}{
Optional Boolean, set to True by default and ignored if 'threshold' is given. When TRUE, the maximal provided FISH distance will be used for the distance range threshold. Otherwise, the threshold will be estimated by the maximal FISH distance that present a small enough deviation (< delta) from the model. 
}
\item{delta}{
Optional numeric, set to 0.05 by default and ignored if 'threshold' is given or if 'useMax' is set to TRUE. Defines the acceptable deviation from the model, when the distance range threshold is estimated from the fit (see details).
}
\item{buffer}{
Optional numeric, set to 1.0 by default and ignored if 'useMax' is set to FALSE. Defines a constant that is added to the threshold value when 'useMax' is set to TRUE. 
}
}

\details{
We use a power law model to relate a set of FISH distances, D, and a matching set of contact frequencies, C:
\eqn{C} ~ \eqn{\beta}D^\eqn{\alpha} \cr 
Taking the log of this equation gives a linear dependency:
\eqn{log(C)}~ \eqn{log(\beta)} + \eqn{\alpha}\eqn{log(D)} \cr                                                          	
Here, we consider only a subset of distances for solving the latter equation and estimate alpha and beta with a linear regression. 
The threshold t, defining the range limit of Hi-C (a distance above which Hi-C frequencies are no longer informative) could be set to the maximal distance in D, or estimated more restrictively from the fit:\cr 
t = maxD\{ | e^( (log(C)-log(\eqn{\beta)})/\eqn{\alpha} ) - D |< \eqn{\delta} \} }                                                                   

\value{
A list with the following objects:
\item{calib}{a list defining the calibration, with the following objects:
f - the calibration function (the power law model), and params - a list of parameters for f (the parameters of the model and the threshold).
}
\item{fit}{the return value of \code{\link{lm}}, used to solve the linear regression}
}
\references{
Y. Shavit, F.K. Hamey, P. Lio', FisHiCal: an R package for iterative FISH-based calibration of Hi-C data, 2014 (submitted).
}
\author{
Yoli Shavit
}
\seealso{
\code{\link{prepareData}}\cr
\code{\link{calibrate}}\cr
\code{\link{conf}}\cr
\code{\link{hic}}\cr
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
 fit = res$fit
 alpha = calib$params[[1]]
 beta = calib$params[[2]]
 threshold = calib$params[[3]]
 
 # plot
 plot(match$frequencies ~ match$distances, xlab = "distances", 
 ylab = "frequencies")
 
 lines((exp(beta)*match$distances^alpha)~match$distances, 
 col = "red")
 
 plot(log(match$frequencies) ~ log(match$distances), 
 xlab = "log(distances)", ylab = "log(frequencies)")
 abline(fit, col = "red")
 
 # plot the estimated threshold 
 abline(h = beta + log(threshold)*alpha, lty = 3)
             
}


