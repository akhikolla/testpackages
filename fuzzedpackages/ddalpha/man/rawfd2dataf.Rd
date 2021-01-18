\name{rawfd2dataf}
\alias{rawfd2dataf}
\title{Transform Raw Functional Data to a \code{dataf} Object}
\usage{
rawfd2dataf(X, range)
}
\arguments{
\item{X}{Either a matrix of size \code{n*d}, or an array of dimension \code{n*d*k} of functional values. Here \code{n}
stands for the number of functions, \code{d} is the number of equi-distant points in the domain where the functional
values are evaluated, and if applicable, \code{k} is the dimensionality of the (vector-valued) functional data.}

\item{range}{A vector of size two that represents the endpoints of the common domain of all functions \code{X}.}
}
\value{
A (possibly multivariate) \code{dataf} object corresponding to the functional data \code{X} evaluated at an
equi-distant grid of points.
}
\description{
Constructs a (possibly multivariate) functional data object given by an array of its functional values
evaluated at an equi-distant grid of points, and transforms it into a \code{dataf} object more suitable 
for work in the \code{ddalpha} package.
}
\examples{
## transform a matrix into a functional data set
n = 5
d = 21
X = matrix(rnorm(n*d),ncol=d)
rawfd2dataf(X,range=c(0,1))

## transform an array into a multivariate functional data set
k = 3
X = array(rnorm(n*d*k),dim=c(n,d,k))
rawfd2dataf(X,range=c(-1,1))

}
\seealso{
\code{\link{dataf2rawfd}}

\code{\link{depthf.fd1}}

\code{\link{depthf.fd2}}
}
\author{
Stanislav Nagy, \email{nagy at karlin.mff.cuni.cz}
}
\keyword{functional}
