\name{VARXpred}
\alias{VARXpred}
\title{VARX Model Prediction
}
\description{Computes point forecasts of a VARX model. The values of 
exogenous variables must be given.
}
\usage{
VARXpred(m1, newxt = NULL, hstep = 1, orig = 0)
}
\arguments{
  \item{m1}{An output object of VARX or refVARX command
}
  \item{newxt}{The data matrix of exogenous variables needed in forecasts.
}
  \item{hstep}{Forecast horizon
}
  \item{orig}{Forecast origin. Default is 0, meaning the last data point.
}
}
\details{Uses the provided exogenous variables and the model to 
compute forecasts
}
\value{Point forecasts and their standard errors
}
\references{Tsay (2014, Chapter 6). Multivariate Time Series Analysis with 
R and Financial Applications. John Wiley. Hoboken, NJ.
}
\author{Ruey S. Tsay
}




