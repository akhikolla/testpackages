\name{REGtspred}
\alias{REGtspred}
\title{Prediction of a fitted regression model with time series errors
}
\description{Perform prediction of a REGts model
}
\usage{
REGtspred(model,newxt,h=1,orig=0)
}
\arguments{
  \item{model}{An output of the REGts command for a vector time series with exogenous variables
}
  \item{newxt}{The new data matrix of the exogenous variables. It must be of the same dimension as the 
original exogenous variables and of length at least h (the forecast horizon).
}
  \item{orig}{The forecast origin. The default is zero indicating that the origin is the 
last observation.
}
 \item{h}{The forecast horizon. For a given h, it computes 1-step to h-step ahead forecasts. Default is 1.
 }
}
\details{Perform prediction of a fitted REGts model
}
\value{
\item{pred }{Forecasts}
\item{se.err }{Standard errors of forecasts}
\item{rmse }{Root mean squares of forecast errors}
\item{rmse }{Root mean squared forecast errors}
\item{orig}{Return the forecast origin}
}
\references{Tsay (2014). Multivariate Time Series Analysis with R and Financial 
Applications. John Wiley. Hoboken, NJ.
}
\author{Ruey S. Tsay
}