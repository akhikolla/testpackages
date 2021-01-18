\name{VARMAcov}
\alias{VARMAcov}
\title{Autocovariance Matrices of a VARMA Model
}
\description{Uses psi-weights to compute the autocovariance matrices of 
a VARMA model
}
\usage{
VARMAcov(Phi = NULL, Theta = NULL, Sigma = NULL, lag = 12, trun = 120)
}
\arguments{
  \item{Phi}{A k-by-kp matrix consisting of VAR coefficient matrices, 
Phi = [Phi1, Phi2, ..., Phip].
}
  \item{Theta}{A k-by-kq matrix consisting of VMA coefficient matrices,
Theta = [Theta1, Theta2, ..., Thetaq]
}
  \item{Sigma}{Covariance matrix of the innovations (k-by-k).
}
  \item{lag}{Number of cross-covariance matrices to be computed. 
Default is 12.
}
  \item{trun}{The lags of pis-weights used in calculation. Default is 120.
}
}
\details{Use psi-weight matrices to compute approximate autocovariance 
matrices of a VARMA model.
}
\value{
\item{autocov}{Autocovariance matrices}
\item{ccm}{Auto correlation matrices}
}
\references{Tsay (2014, Chapter 3). Multivariate Time Series Analysis with R and 
Financial Applications. John Wiley. Hoboken, NJ.
}
\author{Ruey S. Tsay
}
\examples{
Phi=matrix(c(0.2,-0.6,0.3,1.1),2,2)
Sig=matrix(c(4,1,1,1),2,2)
VARMAcov(Phi=Phi,Sigma=Sig)
}




