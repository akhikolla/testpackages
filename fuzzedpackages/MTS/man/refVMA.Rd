\name{refVMA}
\alias{refVMA}
\title{Refining VMA Models
}
\description{Refines a fitted VMA model by setting insignificant parameters to zero
}
\usage{
refVMA(model, thres = 1)
}
\arguments{
  \item{model}{An output object from the command VMA or the refVMA command
}
  \item{thres}{A threshold for individual t-ratio of parameter estimate. 
Default is 1.
}
}
\details{The program simultaneously sets all estimates with t-ratios 
less than the threshold (in modulus) to zero.
}
\value{The same as those of the command VMA.
}
\references{Tsay (2014, Chapter 3). Multivariate Time Series Analysis with R and 
Financial Applications. John Wiley. Hoboken, NJ.
}
\author{Ruey S. Tsay
}
\seealso{VMA
}




