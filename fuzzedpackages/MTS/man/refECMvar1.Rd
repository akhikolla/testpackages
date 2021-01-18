\name{refECMvar1}
\alias{refECMvar1}
\title{Refining ECM for a VAR process
}
\description{Performs constrained least squares estimation of 
a ECM VAR model with known co-integrated processes
}
\usage{
refECMvar1(m1, thres = 1)
}
\arguments{
  \item{m1}{An object of the ECMvar1 command or the refECMvar1 command
}
  \item{thres}{Threshold for individual t-ratio. Default is 1.
}
}
\details{Setting all estimates with t-ration less than the threshold, 
in absoluate value, to zero simultaneously.
}
\value{Constrained estimation results of an ECM VAR model
}
\references{Tsay (2014, Chapter 5). Multivariate Time Series Analysis 
with R and Financial Applications. John Wiley. Hoboken, NJ.
}
\author{Ruey S. Tsay
}
\seealso{ECMvar1, refECMvar
}



