\name{VARMAirf}
\alias{VARMAirf}
\title{Impulse Response Functions of a VARMA Model
}
\description{Compute and plot the impulse response function of a given VARMA model
}
\usage{
VARMAirf(Phi = NULL, Theta = NULL, Sigma = NULL, lag = 12, orth = TRUE)
}
\arguments{
  \item{Phi}{A k-by-kp matrix of VAR coefficients in the form 
Phi=[Phi1, Phi2, ..., Phip].
}
  \item{Theta}{A k-by-kq matrix of VMA coefficients in the form 
Theta=[Theta1, Theta2, ..., Thetaq]
}
  \item{Sigma}{Covariance matrix (k-by-k) of the innovations.
}
  \item{lag}{Number of lags of impulse response functions to be computed
}
  \item{orth}{A logical switch to use orthogonal innovations.
Deafult is to perform orthogonalization of the innovations.
}
}
\value{
\item{psi}{The Psi-weight matrices}
\item{irf}{Impulse response functions}
}
\references{Tsay (2014, Chapter 3). Multivariate Time Series Analysis with 
R and Financial Applications. John Wiley. Hoboken, NJ.
}
\author{Ruey S. Tsay 
}
\seealso{VARMApsi command
}
\examples{
p1=matrix(c(0.2,-0.6,0.3,1.1),2,2)
th1=matrix(c(-0.5,0.2,0.0,-0.6),2,2)
Sig=matrix(c(4,1,1,1),2,2)
m1=VARMAirf(Phi=p1,Theta=th1,Sigma=Sig)
}




