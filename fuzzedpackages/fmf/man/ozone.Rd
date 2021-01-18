\name{ozone}
\alias{ozone}
\docType{data}
\title{
Ozone Level Detection Data Set
}
\description{
Forecasting skewed biased stochastic ozone days: analyses, solutions and beyond, Knowledge and Information Systems, Vol. 14, No. 3, 2008.
Discusses details about the dataset, its use as well as various experiments (both cross-validation and streaming) using many state-of-the-art methods.

A shorter version of the paper (does not contain some detailed experiments as the journal paper above) is in:
Forecasting Skewed Biased Stochastic Ozone Days: Analyses and Solutions. ICDM 2006: 753-764 

}
\usage{data(ozone)}
\format{
  A data frame with 2536 Instances and 73 attributes (including the class attribute, "Class": ozone day, normal day)
}
\details{
The following are specifications for several most important attributes 
that are highly valued by Texas Commission on Environmental Quality (TCEQ). 
More details can be found in the two relevant papers.
 
-- O 3 - Local ozone peak prediction
-- Upwind - Upwind ozone background level
-- EmFactor - Precursor emissions related factor
-- Tmax - Maximum temperature in degrees F
-- Tb - Base temperature where net ozone production begins (50 F)
-- SRd - Solar radiation total for the day
-- WSa - Wind speed near sunrise (using 09-12 UTC forecast mode)
-- WSp - Wind speed mid-day (using 15-21 UTC forecast mode) 

}
\source{
Kun Zhang
zhang.kun05 '@' gmail.com
Department of Computer Science, 
Xavier University of Lousiana

Wei Fan
wei.fan '@' gmail.com
IBM T.J.Watson Research

XiaoJing Yuan
xyuan '@' uh.edu
Engineering Technology Department, 
College of Technology, University of Houston 
}
\references{
https://archive.ics.uci.edu/ml/datasets/Ozone+Level+Detection
}
\examples{

data(ozone)

X = ozone[, -1]

y = ozone[, 1]
}
\keyword{datasets}
