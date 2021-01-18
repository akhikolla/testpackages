\name{FHCSpsdo}
\alias{FHCSpsdo}
\title{Computes the univariate MCD estimator of scatter}
\usage{
  FHCSpsdo(z0,h=NULL,seed=1,q=NULL,ndir=1000)
}
\arguments{
  \item{z0}{Either a data matrix or the result of a call to \code{FHCSkernelEVD}.}
	\item{h}{Number of observation used to compute the univairate outlyingness. Defaults to \code{[(n+q+1)/2]+1}.}
 	\item{seed}{Seed used to initialize the RNG. Defaults to 1.}
	 \item{q}{Number of components. Defaults to \code{ncol(z0)}.}
	 \item{ndir}{Number of projection used to compute the PP outlyngness.}
}
\value{
   A  list with components:
	 	\item{rawDist:}{Outlyingness index of the data on the raw q-dimensonal subset that initialized H*.}
		\item{best:}{the indexes of the members of the H+, the FastHSC subset after the C-steps.} 
		\item{center:}{the p-vector of column means of the observations with indexes in \code{best}.}
		\item{loadings:}{the (rank q) loadings matrix of the observations with indexes in \code{best}.}
		\item{eigenvalues:}{the first \code{min(q)} eigenvalues of the observations with indexes in \code{best}.}
}
\description{
 Pseudo Stahel Donoho Outlyingness based estimates of PCA.
}
\examples{
n<-50
p<-10
x<-matrix(rnorm(n*p),nc=p)
FHCSpsdo(x)
}
\author{
Vakili Kaveh.
}
\references{
	 Rousseeuw, P. J. (1984), Least Median of Squares Regression,
Journal of the American Statistical Association,79,871--880.
}
\keyword{multivariate}
\keyword{robust}
