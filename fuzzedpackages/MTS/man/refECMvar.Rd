\name{refECMvar}
\alias{refECMvar}
\title{Refining Error-Correction Model for VAR series
}
\description{Refining an estimated ECM VAR model by setting insignificant 
estimates to zero
}
\usage{
refECMvar(m1, thres = 1)
}
\arguments{
  \item{m1}{An object of the ECMvar command or the refECMvar command
}
  \item{thres}{Threshold for individual t-ratio. The default is 1.
}
}
\details{Set simultaneously all estimates with t-ratio less than 
the threshold to zero (in modulus).
}
\value{Constrained estimation results of a ECM VAR model
}
\references{Tsay (2014, Chapter 5). Multivariate Time Series Analysis with 
R and Financial Applications. John Wiley. Hoboken, NJ.
}
\author{Ruey S. Tsay
}




