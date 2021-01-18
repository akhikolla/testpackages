\name{tfm2}
\alias{tfm2}
\title{Transfer Function Model with Two Input Variables
}
\description{Estimation of a general transfer function model with two input variables. The model can 
handle one output and up-to 2 input variables. The time series noise can assume multiplicative 
seasonal ARMA models.
}
\usage{
tfm2(y,x,x2=NULL,ct=NULL,wt=NULL,orderN=c(1,0,0),orderS=c(0,0,0),
sea=12,order1=c(0,1,0),order2=c(0,-1,0))
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
 \item{sea}{Seasonality, default is 12 for monthly data
}
  \item{order1}{Order (r,s,b) of the transfer function model of the first input variable, where 
r and s are the degrees of denominator and numerator polynomials and 
b is the delay
}
 \item{order2}{Order (r2,s2,b2) of the transfer function model of the second input variable, where 
2r and s2 are the degrees of denominator and numerator polynomials and 
b2 is the delay
}
}
\details{Perform estimation of a general transfer function model with two input variables
}
\value{
\item{estimate }{Coefficient estimates}
\item{sigma2 }{Residual variance sigma-square}
\item{residuals }{Residual series}
\item{varcoef }{Variance of the estimates}
\item{Nt }{The disturbance series}
\item{rAR}{Regular AR coefficients}
\item{rMA}{Regular MA coefficients}
\item{sAR}{Seasonal AR coefficients}
\item{sMA}{Seasonal MA coefficients}
\item{omega}{Numerator coefficients of the first transfer function}
\item{delta}{Denominator coefficients of the first transfer function}
\item{omega2}{Numerator coefficients of the 2nd transfer function}
\item{delta2}{Denominator coefficients of the 2nd transfer function}
}
\references{Box, G. E. P., Jenkins, G. M., and Reinsel, G. C. (1994). Time Series Analysis:
Forecasting and Control, 3rd edition, Prentice Hall, Englewood Cliffs, NJ.
}
\author{Ruey S. Tsay
}
\seealso{tfm, tfm1
}



