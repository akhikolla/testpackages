\name{VARorderI}
\alias{VARorderI}
\title{VAR order specification I
}
\description{This program is similar to VARorder, but it uses 
observations from t=p+1 to T to compute the information criteria for 
a given VAR(p) model.
}
\usage{
VARorderI(x, maxp = 13, output = T)
}
\arguments{
  \item{x}{A T-by-k data matrix of vector time series
}
  \item{maxp}{The maximum VAR order entertained
}
  \item{output}{A logical switch to control output
}
}
\details{For a given VAR(p) model, the program uses observations 
from t=p+1 to T to compute the information criteria. Therefore, 
different numbers of data points are used to estimate different 
VAR models.
}
\value{
\item{aic}{Akaike information criterion}
\item{aicor}{Order selected by AIC}
\item{bic}{Bayesian information criterion}
\item{bicor}{Order selected by BIC}
\item{hq}{Hannan and Quinn information criterion}
\item{hqor}{Order selected by hq}
\item{Mstat}{Step-wise Chi-square statistics}
\item{Mpv}{p-values of the M-statistics}
}
\references{Tsay (2014)
}
\author{Ruey S Tsay
}
\seealso{VARorder
}




