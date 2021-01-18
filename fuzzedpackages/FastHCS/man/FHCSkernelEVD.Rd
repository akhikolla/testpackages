\name{FHCSkernelEVD}
\alias{FHCSkernelEVD}
\title{Carries out the kernelEVD algorithm for data reduction}
\usage{
  FHCSkernelEVD(x,best=NULL,q=NULL)
}
\arguments{
  \item{x}{A data matrix.}
 \item{best}{An optional subset of \code{1:n}.}
  \item{q}{Desired rank of the SVD decomposition. Optional.}

}
\value{
  A reduced data set with full rank.
}
\description{
This step reduces the data space to the affine subspace spanned by the
\code{n} observations. 
}
\examples{
n<-50
p<-200
x<-matrix(rnorm(n*p),nc=p)
W<-FHCSkernelEVD(x)
}
\author{
Small modification of the code from the \code{classPC} from \code{rrcov}. 
}
\references{
	 Wu, W., Massart, D. L., and de Jong, S. (1997), 'The Kernel PCA Algorithms
        for Wide Data. Part I: Theory and Algorithms,' Chemometrics and Intelligent
        Laboratory Systems,36,165--172
}
\keyword{multivariate}
\keyword{robust}
