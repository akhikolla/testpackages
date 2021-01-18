\name{comVol}
\alias{comVol}
\title{Common Volatility
}
\description{Compute the principal volatility components based on the 
residuals of a VAR(p) model.
}
\usage{
comVol(rtn, m = 10, p = 1, stand = FALSE)
}
\arguments{
  \item{rtn}{A T-by-k data matrix of k-dimensional asset returns
}
  \item{m}{The number of lags used to compute generalized cross-Kurtosis 
matrix
}
  \item{p}{VAR order for the mean equation
}
  \item{stand}{A logical switch to standardize the returns
}
}
\details{Perform a VAR(p) fit, if any. Then, use the residual series to 
perform principal volatility component analysis. The ARCH test statistics 
are also computed for the sample principal components
}
\value{
\item{residuals}{The residuals of a VAR(p) fit}
\item{values}{Eigenvalues of the principal volatility component analysis}
\item{vectors}{Eigenvectors of the principal volatility component analysis}
\item{M}{The transformation matrix}
}
\references{Tsay (2014, Chapter 7)
}
\author{Ruey S. Tsay and Y.B. Hu
}
\examples{
data("mts-examples",package="MTS")
zt=diffM(log(qgdp[,3:5]))
m1=comVol(zt,p=2)
names(m1)
}




