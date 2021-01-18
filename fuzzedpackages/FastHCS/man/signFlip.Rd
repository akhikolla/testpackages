\name{signFlip}
\alias{signFlip}
\title{Carries out the signflip adjustment of a loadings matrix}
\usage{
  signFlip(loadings)
}
\arguments{
  \item{loadings}{A matrix of loadings.}

}
\value{
  An (eventually sign flipped) loadings matrix.
}
\description{
This function solves the sign indeterminancy of the loadings by 
setting the maximum element in a singular vector to be positive. 
}
\examples{
x<-diag(10)
x[1,1]<--2
W<-signFlip(x)
W[1,1]
}
\author{
Kaveh Vakili
}
\keyword{multivariate}
\keyword{robust}
