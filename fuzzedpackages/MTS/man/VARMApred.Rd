\name{VARMApred}
\alias{VARMApred}
\title{VARMA Prediction
}
\description{Compute forecasts and their associate forecast error covariances 
of a VARMA model
}
\usage{
VARMApred(model, h = 1, orig = 0)
}
\arguments{
  \item{model}{A fitted VARMA model
}
  \item{h}{Number of steps of forecasts, i.e., forecast horizon.
}
  \item{orig}{Forecast origin. Default is the end of the sample.
}
}
\value{
\item{pred}{Predictions}
\item{se.err}{Standard errors of forecasts}
\item{orig}{Forecast origin}
}
\references{Tsay (2014, Chapter 3). Multivariate Time Series 
Analysis with R and Financial Applications. John Wiley. 
Hoboken, NJ.
}
\author{Ruey S. Tsay
}




