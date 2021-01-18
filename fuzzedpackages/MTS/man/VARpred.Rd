\name{VARpred}
\alias{VARpred}
\title{VAR Prediction
}
\description{Computes the forecasts of a VAR model, the associated standard errors of forecasts and 
the mean squared errors of forecasts
}
\usage{
VARpred(model, h = 1, orig = 0, Out.level = FALSE, output = TRUE)
}
\arguments{
  \item{model}{An output object of a VAR or refVAR command
}
  \item{h}{Forecast horizon, a positive integer
}
  \item{orig}{Forecast origin. Default is zero meaning the forecast origin 
is the last data point
}
  \item{Out.level}{Boolean control for details of output
}
  \item{output}{Boolean control for printing forecast results
}
}
\details{Computes point forecasts and the associated variances of forecast errors
}
\value{
\item{pred}{Point predictions}
\item{se.err}{Standard errors of the predictions}
\item{mse}{Mean-square errors of the predictions}
}
\references{Tsay (2014, Chapter 2). Multivariate Time Series Analysis with R and 
Financial Applications. John Wiley. Hoboken, NJ.
}
\author{Ruey S. Tsay
}
\examples{
data("mts-examples",package="MTS")
gdp=log(qgdp[,3:5])
zt=diffM(gdp)
m1=VAR(zt,p=2)
VARpred(m1,4)
}




