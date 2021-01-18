\name{msqrt}
\alias{msqrt}
\title{Square Root Matrix
}
\description{Compute the symmetric square root of a positive definite 
matrix
}
\usage{
msqrt(M)
}
\arguments{
  \item{M}{A positive definite matrix
}
}
\details{Use spectral decomposition to compute the square root of 
a positive definite matrix
}
\value{
\item{mtxsqrt}{The square root matrix}
\item{invsqrt}{The inverse of the square root matrix}
}
\author{Ruey S. Tsay
}
\note{This command is used in some of the MTS functions.
}
\examples{
m=matrix(c(1,0.2,0.2,1),2,2)
m1=msqrt(m)
names(m1)
}




