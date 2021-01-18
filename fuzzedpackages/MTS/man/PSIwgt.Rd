\name{PSIwgt}
\alias{PSIwgt}
\title{Psi Wights Matrices
}
\description{Computes the psi-weight matrices of a VARMA model
}
\usage{
PSIwgt(Phi = NULL, Theta = NULL, lag = 12, plot = TRUE, output = FALSE)
}
\arguments{
  \item{Phi}{A k-by-kp matrix of VAR coefficient matrix. 
Phi=[Phi1, Phi1, ..., Phip]
}
  \item{Theta}{A k-by-kq matrix of VMA coefficient matrix. 
Theta=[Theta1, Theta2, ..., Thetaq]
}
  \item{lag}{The number of psi-weight matrices to be computed. 
Deafult is 12.
}
  \item{plot}{A logical switch to control plotting of the psi-weights.
}
  \item{output}{A logical switch to control the output.
}
}
\value{
\item{psi.weight}{Psi-weight matrices}
\item{irf}{Impulse response cofficient matrices}
}
\references{Tsay (2014, Chapter 3). Multivariate Time Series Analysis with R and 
Financial Applications. John Wiley. Hoboken, NJ.
}
\author{Ruey S. Tsay
}
\examples{
phi=matrix(c(0.2,-0.6,0.3,1.1),2,2)
theta=matrix(c(-0.5,0.2,0.0,-0.6),2,2)
m1=PSIwgt(Phi=phi,Theta=theta)
}




