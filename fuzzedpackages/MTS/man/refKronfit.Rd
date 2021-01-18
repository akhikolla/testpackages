\name{refKronfit}
\alias{refKronfit}
\title{Refining VARMA Estimation via Kronecker Index Approach
}
\description{This program performs model simplification of a fitted 
VARMA model via the Kronecker index approach
}
\usage{
refKronfit(model, thres = 1)
}
\arguments{
  \item{model}{The name of a model from the command Kronfit or refKronfit
}
  \item{thres}{A threshold for t-ratio of individual parameter estimate.
The default is 1.
}
}
\details{For a given threshold, the program sets a parameter to zero if 
its t-ratio (in modulus) is less than the threshold.
}
\value{Same as those of the command Kronfit.
}
\references{Tsay (2014, Chapter 4). Multivariate Time Series Analysis with R 
and Financial Applications. John Wiley. Hoboken, NJ.
}
\author{Ruey S. Tsay
}
\seealso{Kronfit
}




