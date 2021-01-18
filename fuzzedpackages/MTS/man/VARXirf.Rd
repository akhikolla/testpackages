\name{VARXirf}
\alias{VARXirf}
\title{Impluse response function of a fitted VARX model
}
\description{Compute the impulse response functions and cumulative impulse response functions of 
a fitted VARX model
}
\usage{
VARXirf(model,lag=12,orth=TRUE)
}
\arguments{
  \item{model}{An output of the VARX (or refVARX) command for a vector time series with exogeneous variables
}
  \item{lag}{The number of lags of the impulse response function to be computed. Default is 12.
}
  \item{orth}{The control variable for using orthogonal innovations. This command applies 
to the impulse response functions of the VAR part only.
}
}
\details{Compute the impulse response functions and cumulative impulse response functions 
of a fitted VARX model. The impulse response function of the exogeneous variables are 
also given. The plots of impulse response functions are shown.
}
\value{
\item{irf }{Impulse response functions of the VAR part, original innovations used}
\item{orthirf }{Impulse response functions of the VAR part using orthogonal innovations}
\item{irfX}{Impulse response function of the exogenous variables}
}
\references{Tsay (2014). Multivariate Time Series Analysis with R and Financial 
Applications. John Wiley. Hoboken, NJ.
}
\author{Ruey S. Tsay
}
