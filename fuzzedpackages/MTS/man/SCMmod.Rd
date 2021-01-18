\name{SCMmod}
\alias{SCMmod}
\title{Scalar Component Model specification
}
\description{For a given set of SCMs and locator of transformation matrix, 
the program specifies a VARMA model via SCM approach for estimation
}
\usage{
SCMmod(order, Ivor, output)
}
\arguments{
  \item{order}{A k-by-2 matrix of the orders of SCM
}
  \item{Ivor}{A k-dimensional vector indicating the location of "1" 
for each component in the transformation matrix.
}
  \item{output}{A logical switch to control output.
}
}
\details{The command specified estimable parameters for a VARMA model via 
the SCM components. In the output, "2" denotes estimation, "1" denotes fixing 
the value to 1, and "0" means fixing the parameter to zero.
}
\value{
\item{Tmtx}{Specification of the transformation matrix T}
\item{ARpar}{Specification of the VAR parameters}
\item{MApar}{Specification of the VMA parameters}
}
\references{Tsay (2014, Chapter 4). Multivariate Time Series Analysis with R and 
Financial Applications. John Wiley. Hoboken, NJ.
}
\author{Ruey S. Tsay
}
\examples{
ord=matrix(c(0,1,1,0,0,1),3,2)
Ivor=c(3,1,2)
m1=SCMmod(ord,Ivor,TRUE)
}




