\name{RealVAMS-package}
\alias{RealVAMS-package}
\docType{package}
\title{
Multivariate VAM Fitting
}
\description{
Fits a multivariate value-added model (VAM), see Broatch, Green, and Karl (2018) \cr <doi:10.32614/RJ-2018-033>, and Broatch and Lohr (2012) \cr <doi:10.3102/1076998610396900>, with normally distributed test scores and a binary outcome indicator. A pseudo-likelihood approach, Wolfinger and O'Connell (1993) \cr <doi:10.1080/00949659308811554>, is used for the estimation of this joint generalized linear mixed model.  The inner loop of the pseudo-likelihood routine (estimation of a linear mixed model) occurs in the framework of the EM algorithm presented by \cr Karl, Yang, and Lohr (2013) <DOI:10.1016/j.csda.2012.10.004>. This material is based upon work supported by the National Science Foundation under grants DRL-1336027 and DRL-1336265.  
}
\details{
\tabular{ll}{
Package: \tab RealVAMS\cr
Type: \tab Package\cr
Version: \tab 0.4-3\cr
Date: \tab 2019-04-22\cr
License: \tab GPL-2\cr
}
}
\author{
Authors: Andrew T. Karl, Jennifer Broatch, and Jennifer Green

Maintainer: Andrew Karl <akarl@asu.edu>
}
\references{
Broatch, J. and Lohr, S. (2012) <DOI:10.3102/1076998610396900> Multidimensional Assessment of Value Added by Teachers to Real-World Outcomes. \emph{Journal of Educational and Behavioral Statistics} \bold{37}, 256--277.

Broatch, J., Green, J., Karl, A. (2018) <DOI:10.32614/RJ-2018-033> RealVAMS: An R Package for Fitting a Multivariate Value-added Model (VAM). \emph{The R Journal} 10/1, 22--30.

Karl, A., Yang, Y. and Lohr, S. (2013) <DOI:10.1016/j.csda.2012.10.004> Efficient Maximum Likelihood Estimation of Multiple Membership Linear Mixed Models, with an Application to Educational Value-Added Assessments. \emph{Computational Statistics & Data Analysis} \bold{59}, 13--27.  

Karl, A., Yang, Y. and Lohr, S. (2013) <DOI:10.3102/1076998613494819> A Correlated Random Effects Model for Nonignorable Missing Data in Value-Added Assessment of Teacher Effects. \emph{Journal of Educational and Behavioral Statistics} \bold{38}, 577--603.

Karl, A., Yang, Y. and Lohr, S. (2014) <DOI:10.1016/j.csda.2013.11.019> Computation of Maximum Likelihood Estimates for Multiresponse Generalized Linear Mixed Models with Non-nested, Correlated Random Effects. \emph{Computational Statistics & Data Analysis} \bold{73}, 146--162.

Lockwood, J., McCaffrey, D., Mariano, L., Setodji, C. (2007) <DOI:10.3102/1076998606298039> Bayesian Methods for Scalable Multivariate Value-Added Assessment. \emph{Journal of Educational and Behavioral Statistics} \bold{32}, 125--150.

Wolfinger, R. (1993) <DOI:10.1080/00949659308811554> Generalized linear mixed models a pseudo-likelihood approach. \emph{Journal of Statistical Computation and Simulation} \bold{48} 233--243.
}

\keyword{ package }

\examples{
data(example.score.data)
data(example.outcome.data)
#The next line exists to show that the function can run and that the package
#installed correctly. This is a CRAN requirement to ensure that the package
#works in future version of R
RealVAMS(example.score.data,example.outcome.data,max.PQL.it=1,max.iter.EM=2,
var.parm.hessian=FALSE)

\donttest{
res<-RealVAMS(example.score.data,example.outcome.data)}
}
