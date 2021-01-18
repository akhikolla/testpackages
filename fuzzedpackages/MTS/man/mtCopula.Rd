\name{mtCopula}
\alias{mtCopula}
\title{Multivariate t-Copula Volatility Model
}
\description{Fits a t-copula to a k-dimensional standardized return series.
The correlation matrices are parameterized by angles and the angles evolve 
over time via a DCC-type equation.
}
\usage{
mtCopula(rt, g1, g2, grp = NULL, th0 = NULL, m = 0, 
				 include.th0 = TRUE, ub=c(0.95,0.049999))
}
\arguments{
  \item{rt}{A T-by-k data matrix of k standardized time series (after 
univariate volatility modeling)
}
  \item{g1}{lamda1 parameter, nonnegative and less than 1
}
  \item{g2}{lambda2 parameter, nonnegative and satisfying lambda1+lambda2 < 1.
}
  \item{grp}{a vector to indicate the number of assets divided into groups.
Default means each individual asset forms a group.
}
  \item{th0}{initial estimate of theta0
}
  \item{m}{number of lags used to estimate the local theta-angles
}
  \item{include.th0}{A logical switch to include theta0 in estimation. 
Default is to include.
}

  \item{ub}{Upper bound of parameters
}
}
\value{
\item{estimates}{Parameter estimates}
\item{Hessian}{Hessian matrix}
\item{rho.t}{Cross-correlation matrices}
\item{theta.t}{Time-varying angel matrices}
}
\references{Tsay (2014, Chapter 7). Multivariate Time Series Analysis 
with R and Financial Applications. John Wiley. Hoboken, NJ.
}
\author{Ruey S. Tsay
}




