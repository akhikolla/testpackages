\name{MTSplot}
\alias{MTSplot}
\title{Multivariate Time Series Plot
}
\description{Provides time plots of a vector time series
}
\usage{
MTSplot(data, caltime = NULL)
}
\arguments{
  \item{data}{data matrix of a vector time series 
}
  \item{caltime}{Calendar time. Default is NULL, that is, using 
time index
}
}
\details{Provides time plots of a vector time series. The output frame 
depends on the dimension of the time series
}
\value{Time plots of vector time series
}
\author{Ruey S. Tsay
}
\examples{
xt=matrix(rnorm(1500),500,3)
MTSplot(xt)
}




