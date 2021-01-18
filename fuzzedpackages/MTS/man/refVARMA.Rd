\name{refVARMA}
\alias{refVARMA}
\title{Refining VARMA Estimation
}
\description{Refines a fitted VARMA model by setting insignificant 
estimates to zero
}
\usage{
refVARMA(model, thres = 1.5)
}
\arguments{
  \item{model}{An output object from the command VARMA or the command 
refVARMA
}
  \item{thres}{A threshold value for individual t-ratio of the estimates.
}
}
\details{The program simultaneously sets estimates with t-ratios less 
than the threshold (in modulus) to zero.
}
\value{The same as those of the command VARMA.
}
\references{Tsay (2014, Chapter 3). Multivariate Time Series Analysis with 
R and Financial Applications. John Wiley. Hoboken, NJ.
}
\author{Ruey S. Tsay
}
\seealso{VARMA
}




