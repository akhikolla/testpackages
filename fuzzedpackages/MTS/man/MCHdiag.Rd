\name{MCHdiag}
\alias{MCHdiag}
\title{Multivariate Conditional Heteroscedastic Model Checking
}
\description{Apply four portmanteau test statistics to check the 
validity of a fitted multivariate volatility model
}
\usage{
MCHdiag(at, Sigma.t, m = 10)
}
\arguments{
  \item{at}{A T-by-k matrix of residuals for a k-dimensional asset return 
series
}
  \item{Sigma.t}{The fitted volatility matrices. The dimension is 
T-by-k^2 matrix
}
  \item{m}{The number of lags used in the tests. Default is 10.
}
}
\details{The four test statistics are given in Tsay (2014, Chapter 7)
}
\value{Four test statistics and their p-values
}
\references{Tsay (2014, Chapter 7). Multivariate Time Series Analysis with 
R and Financial Applications. John Wiley. Hoboken, NJ.
}
\author{Ruey S. Tsay
}




