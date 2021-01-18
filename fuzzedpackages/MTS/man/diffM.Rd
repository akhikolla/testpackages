\name{diffM}
\alias{diffM}
\title{Difference of multivariate time series
}
\description{Performs the difference operation of a vector time series
}
\usage{
diffM(zt, d = 1)
}
\arguments{
  \item{zt}{A vector time series (T by k, with sample size T and dimension k)
}
  \item{d}{Order of differencing. Default is d=1. 
}
}
\details{When d = 1, the command is equivalent to 
apply(zt,2,diff)
}
\value{The differenced time series
}
\author{Ruey S Tsay
}
\examples{
data("mts-examples",package="MTS")
zt=log(qgdp[,3:5])
xt=diffM(zt)
}




