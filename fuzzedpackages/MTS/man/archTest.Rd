\name{archTest}
\alias{archTest}
\title{ARCH test for univariate time series
}
\description{
Perform tests to check the conditional heteroscedasticity in a time series. The Ljung-Box statistics of squared series and a rank-based Ljung-Box 
test are used.
}
\usage{
archTest(rt, lag = 10)
}
\arguments{
  \item{rt}{A scalar time series. If rt is a matrix, only the first column is used.
}
  \item{lag}{
The number of lags of ACF used in the Ljung-Box statistics. The default is 10.
}
}
\details{
The Ljung-Box statistics based on the squared series are computed first. 
The rank series of the squared time series is than used to test the 
conditional heteroscedasticity.
}
\value{The Q-statistic and its p-value. Also, the rank-based Q statistic and its p-value.
}
\author{ Ruey Tsay
}
\seealso{MarchTest
}
\examples{
rt=rnorm(200)
archTest(rt)
}




