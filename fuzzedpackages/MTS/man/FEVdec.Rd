\name{FEVdec}
\alias{FEVdec}
\title{Forecast Error Variance Decomposition
}
\description{Computes the forecast error variance decomposition of a 
VARMA model
}
\usage{
FEVdec(Phi, Theta, Sig, lag = 4)
}
\arguments{
  \item{Phi}{VAR coefficient matrices in the form Phi=[Phi1, Phi2, ..., Phip], 
a k-by-kp matrix.
}
  \item{Theta}{VMA coefficient matrices in form form Theta=[Theta1, Theta2, 
..., Thetaq], a k-by-kq matrix.
}
  \item{Sig}{The residual covariance matrix Sigma, a k-by-k positive definite 
matrix.
}
  \item{lag}{The number of lags of forecast errors variance to be computed. 
Default is 4.
}
}
\details{Use the psi-weight matrices to compute the forecast error covariance 
and use Cholesky decomposition to perform the decomposition
}
\value{
\item{irf}{Impulse response matrices}
\item{orthirf}{Orthogonal impulse response matrices}
\item{Omega}{Forecast error variance matrices}
\item{OmegaR}{Forecast error variance decomposition}
}
\references{Tsay (2014, Chapter 3)
}
\author{Ruey S. Tsay
}
\examples{
p1=matrix(c(0.2,-0.6,0.3,1.1),2,2)
theta1=matrix(c(-0.5,0,0,-0.6),2,2)
Sig=matrix(c(3,1,1,1),2,2)
m1=FEVdec(p1,theta1,Sig)
names(m1)
}




