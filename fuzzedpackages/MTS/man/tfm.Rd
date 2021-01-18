\name{tfm}
\alias{tfm}
\title{Transfer Function Model
}
\description{Estimates a transform function model. This program does not allow 
rational transfer function model. It is a special case of tfm1 and tfm2.
}
\usage{
tfm(y, x, b = 0, s = 1, p = 0, q = 0)
}
\arguments{
  \item{y}{Data vector of dependent (output) variable
}
  \item{x}{Data vector of independent variable
}
  \item{b}{deadtime or delay
}
  \item{s}{The order of the transfer function polynomial
}
  \item{p}{AR order of the disturbance
}
  \item{q}{MA order of the disturbance
}
}
\details{The model entertained is y_t = c_0+ v(B)x_t + n_t. 
v(B) = 1- v1*B - ... - vs*B^s, and n_t is an ARMA(p,q) process.
}
\value{
\item{coef }{Coefficient estimates of the transfer function}
\item{se.coef }{Standard errors of the transfer function coefficients}
\item{coef.arma }{Coefficient estimates of ARMA models}
\item{se.arma }{Standard errors of ARMA coefficients}
\item{nt }{The disturbance series}
\item{residuals }{The residual series}
}
\references{Box, G. E. P., Jenkins, G. M., and Reinsel, G. C. (1994). 
Time Series Analysis:
Forecasting and Control, 3rd edition, Prentice Hall, Englewood Cliffs, NJ.
}
\author{Ruey S. Tsay
}




