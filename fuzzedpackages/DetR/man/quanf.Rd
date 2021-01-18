\name{quanf}
\alias{quanf}
\title{Converts alpha values to h-values}
\usage{
  quanf(n,p,alpha)
}
\arguments{
 \item{n}{Number of rows of the data matrix.}
  \item{p}{Number of columns of the data matrix.}
\item{alpha}{Numeric parameter controlling the size of the active subsets,
	 i.e., \code{"h=quanf(alpha,n,p)"}.  Allowed
          values are between 0.5 and 1 and the default is 0.5.}
}
\value{
  An integer number of the size of the starting p-subsets.
}
\description{
  DetLTS selects the subset of size h that minimizes the log-determinant criterion. 
The function \code{quanf} determines the size of h based on the rate of contamination
 the user expects is present in the data.
  This is an internal function not intended to be called 
  by the user.
}
\examples{
quanf(p=3,n=500,alpha=0.5)
}
\author{
Kaveh Vakili
}
\keyword{multivariate}
\keyword{robust}

