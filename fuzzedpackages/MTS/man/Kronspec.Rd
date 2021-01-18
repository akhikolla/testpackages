\name{Kronspec}
\alias{Kronspec}
\title{Kronecler Index Specification
}
\description{For a given set of Kronecker indices, the program 
specifies a VARMA model. It gives details of parameter specification.
}
\usage{
Kronspec(kdx, output = TRUE)
}
\arguments{
  \item{kdx}{A vector of Kronecker indices
}
  \item{output}{A logical switch to control output.  Default is with 
output.
}
}
\value{
\item{PhiID}{Specification of the AR matrix polynomial. 0 denotes zero 
parameter, 1 denotes fixing parameter to 1, and 2 denotes the parameter 
requires estimation}
\item{ThetaID}{Specification of the MA matrix polynomial}
}
\references{Tsay (2014, Chapter 4)
}
\author{Ruey S. Tsay
}
\examples{
kdx=c(2,1,1)
m1=Kronspec(kdx)
names(m1)
}




