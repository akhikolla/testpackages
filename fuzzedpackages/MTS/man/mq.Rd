\name{mq}
\alias{mq}
\title{Multivariate Ljung-Box Q Statistics
}
\description{Computes the multivariate Ljung-Box statistics for 
cross-correlation matrices
}
\usage{
mq(x, lag = 24, adj = 0)
}
\arguments{
  \item{x}{The data matrix of a vector time series or residual series of a 
fitted multivariate model.
}
  \item{lag}{The number of cross-correlation matrices used. Default is 24.
}
  \item{adj}{Adjustment for the degrees of freedom for the Ljung-Box 
statistics. This is used for residual series. Default is zero.
}
}
\details{Computes the multivariate Ljung-Box statistics and their 
p-values. For model checking, the subcommand adj can be used to 
adjust the degrees of freedom of the Chi-square statistics.
}
\value{The multivariate Q-statistics and their p-values. Also, it provides a 
plot of the p-values.
}
\references{Tsay (2014). Multivariate Time Series Analysis with R and 
Financial Applications. John Wiley. Hoboken, NJ.
}
\author{Ruey S. Tsay
}
\examples{
x=matrix(rnorm(1500),500,3)
mq(x)
}




