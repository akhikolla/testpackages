\name{PIwgt}
\alias{PIwgt}
\title{Pi Weight Matrices
}
\description{Compute the Pi-weight matrices of a VARMA model
}
\usage{
PIwgt(Phi = NULL, Theta = NULL, lag = 12, plot = TRUE)
}
\arguments{
  \item{Phi}{A k-by-kp matrix of VAR coefficients in the form 
[Phi1, Phi2, Phi3, ..., Phip]
}
  \item{Theta}{A k-by-kq matrix of VMA coefficients in the form 
[Theta1, Theta2, ..., Thetaq]
}
  \item{lag}{The number of Pi-weight matrices to be computed.
}
  \item{plot}{A logical switch to plot the Pi-weight matrices
}
}
\details{The Pi-weight matrices for a VARMA model is 
Pi(B) = inverse(Theta(B)) times Phi(B). 
}
\value{
\item{pi.weight}{The matrix of Pi-weight coefficient} 
}
\references{Tsay (2014, Chapters 2 and 3). 
Multivariate Time Series Analysis with R and Financial Applications. 
John Wiley. Hoboken, NJ.
}
\author{Ruey S. Tsay
}
\seealso{PSIwgt
}
\examples{
Phi1=matrix(0,2,2); Phi2=matrix(c(0.2,-0.6,0.3,1.1),2,2)
Theta1=diag(c(-0.5,-0.4))
Phi=cbind(Phi1,Phi2)
m1=PIwgt(Phi=Phi,Theta=Theta1)
names(m1)
}




