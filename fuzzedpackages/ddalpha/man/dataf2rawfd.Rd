\name{dataf2rawfd}
\alias{dataf2rawfd}
\title{Transform a \code{dataf} Object to Raw Functional Data}
\usage{
dataf2rawfd(dataf, range = NULL, d = 101)
}
\arguments{
\item{dataf}{Functions to be transformed, represented by a (possibly multivariate) \code{dataf} object of their arguments
and functional values. \code{m} stands for the number of functions. The grid of observation points for the 
functions in \code{dataf} may not be the same.}

\item{range}{The common range of the domain where the functions \code{dataf} are observed.
Vector of length 2 with the left and the right end of the interval. Must contain all arguments given in 
\code{dataf}. If the range is not provided, the smallest interval in which all the arguments from the data functions
are contained is chosen as the domain.}

\item{d}{Grid size to which all the functional data are transformed. All functional observations are 
transformed into vectors of their functional values of length \code{d}
corresponding to equi-spaced points in the domain given by the interval \code{range}. Functional values in these
points are reconstructed using linear interpolation, and extrapolation, see Nagy et al. (2016).}
}
\value{
If the functional data are univariate (scalar-valued), a matrix of size \code{m*d} is given, with each row
corresponding to one function. If the functional data are \code{k}-variate with k>1, an array of size \code{m*d*k}
of the functional values is given.
}
\description{
From a (possibly multivariate) functional data object \code{dataf} constructs an array of the functional values
evaluated at an equi-distant grid of points.
}
\examples{
## transform a matrix into a functional data set and back
n = 5
d = 21
X = matrix(rnorm(n*d),ncol=d)
R = rawfd2dataf(X,range=c(0,1))
R2 = dataf2rawfd(R,range=c(0,1),d=d)
all.equal(X,R2)

## transform a functional dataset into a raw matrix of functional values
dataf = dataf.population()$dataf
dataf2rawfd(dataf,range=c(1950,2015),d=66)

## transform an array into a multivariate functional data set and back
k = 3
X = array(rnorm(n*d*k),dim=c(n,d,k))
R = rawfd2dataf(X,range=c(-1,1))
dataf2rawfd(R,range=c(-1,1),d=50)

}
\seealso{
\code{\link{rawfd2dataf}}

\code{\link{depthf.fd1}}

\code{\link{depthf.fd2}}
}
\author{
Stanislav Nagy, \email{nagy at karlin.mff.cuni.cz}
}
\keyword{functional}
