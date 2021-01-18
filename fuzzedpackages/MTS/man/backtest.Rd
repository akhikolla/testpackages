\name{backtest}
\alias{backtest}
\title{Backtesting of a scalar ARIMA model
}
\description{Perform out-of-sample prediction of a given ARIMA model and compute the summary statistics
}
\usage{
backtest(m1,rt,orig,h,xre=NULL,fixed=NULL,inc.mean=TRUE,reest=1)
}
\arguments{
  \item{m1}{An output of the arima command for scalar time series
}
  \item{rt}{The time series under consideration
}
  \item{orig}{The starting forecast origin. It should be less than the length of the underlying time series
}
 \item{h}{The forecast horizon. For a given h, it computes 1-step to h-step ahead forecasts
 }
  \item{inc.mean}{A logical switch. It is true if mean vector is estimated.
}
  \item{fixed}{A vector of the length of the number of coefficients of the ARIMA model. It is used 
in R for parameter constraint.
}
 \item{xre}{A matrix containing the exogeneous variables used in the ARIMA model
}
 \item{reest}{A control variable used to re-fit the model in prediction. The program will 
re-estimate the model for every new reest observations. The default is 1. That is, re-estimate 
the model with every new data point.
}
}
\details{Perform estimation-prediction-reestimation in the forecasting subsample, and 
to compuate the summary statistics
}
\value{
\item{origion }{Forecast origin}
\item{error }{forecast errors}
\item{forecasts }{forecasts}
\item{rmse }{Root mean squared forecast errors}
\item{mabso}{Mean absolute forecast errors}
\item{reest}{Return the reest value}
}
\references{Tsay (2010). Analysis of Financial Time Series, 3rd. John Wiley. Hoboken, NJ.
}
\author{Ruey S. Tsay
}
