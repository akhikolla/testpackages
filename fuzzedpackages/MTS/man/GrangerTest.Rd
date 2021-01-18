\name{GrangerTest}
\alias{GrangerTest}
\title{Granger Causality Test
}
\description{Performs Granger causality test using a vector autoregressive model
}
\usage{
GrangerTest(X,p=1,include.mean=T,locInput=c(1))
}
\arguments{
  \item{X}{a T-by-p data matrix with T denoting sample size and p the number of 
variables
}
  \item{p}{vector AR order.
}
  \item{include.mean}{Indicator for including a constant in the model. Default is 
TRUE.
}
  \item{locInput}{Locators for the input variables in the data matrix. Default 
is the first column being the input variable. Multiple inputs are allowed.
}
}
\details{Perform VAR(p) and constrained VAR(p) estimations to test the 
Granger causality. It uses likelihood ratio and asymptotic chi-square.
}
\value{
\item{data}{Original data matrix}
\item{cnst}{logical variable to include a constant in the model}
\item{order}{order of VAR model used}
\item{coef}{Coefficient estimates}
\item{constraints}{Implied constraints of Granger causality}
\item{aic, bic, hq}{values of information criteria}
\item{residuals}{residual vector}
\item{secoef}{standard errors of coefficient estimates}
\item{Sigma}{Residual covariance matrix}
\item{Phi}{Matrix of VAR coefficients}
\item{Ph0}{constant vector}
\item{omega}{Estimates of constrained coefficients}
\item{covomega}{covariance matrix of constrained parameters}
\item{locInput}{Locator vector for input variables}
}
\references{Tsay (2014, Chapter 2)
}
\author{Ruey S. Tsay
}





