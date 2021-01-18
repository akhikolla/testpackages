\name{Kronpred}
\alias{Kronpred}
\title{Prediction of a fitted VARMA model via Kronfit, using Kronecker indices
}
\description{Compute forecasts of a fitted VARMA model via the command Kronfit
}
\usage{
Kronpred(model,orig=0,h=1)
}
\arguments{
  \item{model}{A model fitted by the Kronfit command
}
  \item{orig}{Forecast origin. The default is 0, implying that the origin is the 
last observation
}
  \item{h}{Forecast horizon. Default is h=1, 1-step ahead forecast
}
}
\details{For a model, which is the output of the command Kronfit, the command computes 
forecasts of the model starting at the forecast origin. !-step to h-step ahead 
forecasts are computed.
}
\value{
\item{pred}{Forecasts}
 \item{se.err}{Standard errors of the forecasts}
 \item{orig}{Return the forecast origin}
}
\references{Tsay (2014). Multivariate Time Series Analysis with R and Financial Applications,
John Wiley, Hoboken, New Jersey
}
\author{Ruey S. Tsay
}
