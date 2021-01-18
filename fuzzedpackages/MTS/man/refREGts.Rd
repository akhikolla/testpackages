\name{refREGts}
\alias{refREGts}
\title{Refining a Regression Model with Time Series Errors
}
\description{Refines a fitted REGts by setting simultaneously parameters 
with t-ratios less than the threshold (in modulus) to zero
}
\usage{
refREGts(m1, thres = 1)
}
\arguments{
  \item{m1}{An output object from the REGts command or refREGts command
}
  \item{thres}{Threshold value for individual t-ratio. Default is 1.
}
}
\value{The same as those of the command REGts.
}
\references{Tsay (2014, Chapter 6). Multivariate Time Series Analysis with 
R and Financial Applications. John Wiley. Hoboken, NJ.
}
\author{Ruey S. Tsay
}
\seealso{refVAR, refVARMA
}




