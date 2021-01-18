\name{VARpsi}
\alias{VARpsi}
\title{VAR Psi-weights
}
\description{Computes the psi-weight matrices of a VAR model
}
\usage{
VARpsi(Phi, lag = 5)
}
\arguments{
  \item{Phi}{A k-by-kp matrix of VAR coefficients in the form 
Phi=[Phi1, Phi2, ..., Phip]
}
  \item{lag}{Number of psi-weight lags 
}
}
\value{Psi-weights of a VAR model
}
\references{Tsay (2014, Chapter 2). Multivariate Time Series Analysis with R 
and Financial Applications. John Wiley. Hoboken, NJ.
}
\author{Ruey S. Tsay
}
\examples{
p1=matrix(c(0.2,-0.6,0.3,1.1),2,2)
m1=VARpsi(p1,4)
names(m1)
}




