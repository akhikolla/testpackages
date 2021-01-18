\name{sVARMApred}
\alias{sVARMApred}
\title{Prediction of a fitted multiplicative seasonal VARMA model
}
\description{Perform prediction of a seasonal VARMA model
}
\usage{
sVARMApred(model,orig,h=1)
}
\arguments{
  \item{model}{An output of the sVARMA command 
}

  \item{orig}{The forecast origin.
}
 \item{h}{The forecast horizon. For a given h, it computes 1-step to h-step ahead forecasts. Default is 1.
 }
}
\details{Perform prediction of a fitted sVARMA model
}
\value{
\item{data }{The original data matrix}
\item{pred }{Forecasts}
\item{se.err }{Standard errors of forecasts}
\item{orig}{Return the forecast origin}
}
\references{Tsay (2014, chapter 6). Multivariate Time Series Analysis with R and Financial 
Applications. John Wiley. Hoboken, NJ.
}
\author{Ruey S. Tsay
}