\name{MTSdiag}
\alias{MTSdiag}
\title{Multivariate Time Series Diagnostic Checking
}
\description{Performs model checking for a fitted multivariate time series model, including residual cross-correlation matrices, 
multivariate Ljung-Box tests for residuals, and residual plots
}
\usage{
MTSdiag(model, gof = 24, adj = 0, level = F)
}
\arguments{
  \item{model}{A fitted multivariate time series model
}
  \item{gof}{The number of lags of residual cross-correlation matrices used in the tests
}
  \item{adj}{
The adjustment for degrees of freedom of Ljung-Box statistics. 
Typically, the number of fitted coefficients of the model. Default is zero.
}
  \item{level}{Logical switch for printing residual cross-correlation matrices
}
}
\value{Various test statistics, their p-values, and residual plots.
}
\author{Ruey S Tsay
}
\examples{
phi=matrix(c(0.2,-0.6,0.3,1.1),2,2); sigma=diag(2)
m1=VARMAsim(200,arlags=c(1),phi=phi,sigma=sigma)
zt=m1$series
m2=VAR(zt,1,include.mean=FALSE)
MTSdiag(m2)
}




