\name{Btfm2}
\alias{Btfm2}
\title{Back-Test of a Transfer Function Model with Two Input Variables
}
\description{Perform back-test of transfer function model with 2 input variable. 
For a specified tfm2 model and a given forecast origin, the command 
iterated between estimation and 1-step ahead prediction starting at the forecast origin 
until the (T-1)th observation, where T is the sample size.
}
\usage{
Btfm2(y,x,x2=NULL,wt=NULL,ct=NULL,orderN=c(1,0,0),orderS=c(0,0,0),sea=12,
order1=c(0,1,0),order2=c(0,-1,0),orig=(length(y)-1))
}
\arguments{
  \item{y}{Data vector of dependent variable
}
  \item{x}{Data vector of the first input (or independent) variable
}
  \item{x2}{Data vector of the second input variable if any
}
 \item{ct}{Data vector of a given deterministic variable such as time trend, if any
}
 \item{wt}{Data vector of co-integrated series between input and output variables if any
}
  \item{orderN}{Order (p,d,q) of the regular ARMA part of the disturbance component
}
 \item{orderS}{Order (P,D,Q) of the seasonal ARMA part of the disturbance component
}
 \item{sea}{Seasonalityt, default is 12 for monthly data
}
  \item{order1}{Order (r,s,b) of the transfer function model of the first input variable, where 
r and s are the degrees of denominator and numerator polynomials and 
b is the delay
}
 \item{order2}{Order (r2,s2,b2) of the transfer function model of the second input variable, where 
2r and s2 are the degrees of denominator and numerator polynomials and 
b2 is the delay
}
 \item{orig}{Forecast origin with default being T-1, where T is the sample size
}
}
\details{Perform out-of-sample 1-step ahead prediction to evaluate a fitted tfm2 model
}
\value{
\item{ferror }{1-step ahead forecast errors, starting at the given forecast origin}
\item{mse }{out-of-sample mean squared forecast errors}
\item{rmse }{root mean squared forecast errors}
\item{mae }{out-of-sample mean absolute forecast errors}
\item{nobf }{The number of 1-step ahead forecast errors computed}
\item{rAR}{Regular AR coefficients}
}
\references{Box, G. E. P., Jenkins, G. M., and Reinsel, G. C. (1994). Time Series Analysis:
Forecasting and Control, 3rd edition, Prentice Hall, Englewood Cliffs, NJ.
}
\author{Ruey S. Tsay
}
\seealso{tfm2
}



