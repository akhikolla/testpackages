\name{hyperplane}
\alias{hyperplane}
\title{ Hyperplane Passing Through k Points in the k-dimensional Space. }
\description{The function computes a \eqn{(k\!-\!1)}{(k-1)}-dimensional hyperplane passing through \eqn{k} given points in the \eqn{k}-dimensional space.
}
\usage{
hyperplane(X)
}



\arguments{
  \item{X}{ a numeric \eqn{k \times k}{k x k} matrix containing \eqn{k} data points as \strong{rows}.}
}
\details{
  A \eqn{(k\!-\!1)}{(k-1)}-dimensional hyperplane in \eqn{R^k} consists of all points \eqn{x} that satisfy
  \deqn{d^T x + c = 0, }{d' x + c = 0,}
  where \eqn{d} is a \eqn{k}-vector and \eqn{c} is a scalar. The function returns the \eqn{(k\!+\!1)}{(k+1)}-vector \eqn{(d,c)}. 
  It is normalized such that the length of \eqn{d} equals \eqn{(k\!-\!1)!}{(k-1)!} times 
  the \eqn{(k\!-\!1)}{(k-1)}-dimensional volume of the simplex formed by the points on the plane. 
  (If \eqn{k = 3}, this is a triangle.) 
  Hence this function can also easily be used to compute volumes of simplices.
  
  The direction of \eqn{d}, that is, whether it points towards the origin or not, is not fixed. It depends on the
  order of the data points within the matrix \code{X}.   
                                                  
 	If the \eqn{k} points do not uniquely define a \eqn{(k\!-\!1)}{(k-1)}-dimensional hyperplane 
 	(i.e. they lie on a \eqn{(k\!-\!2)}{(k-2)}-dimensional hyperplane), a vector containing zeros is returned. 
}



\value{
  a vector of length \eqn{(k\!+\!1)}{(k-1)} describing the hyperplane, see details above.
}

\author{Daniel Vogel}

\examples{
### ----<< Example 1 >>---- : line in R^2
X <- rbind(c(4,5),c(8,2))
hyperplane(X)
# The line through the the points c(4,5) and c(8,2) is given by 
#    3*x + 4*y - 32 = 0.
# The norm of the first two components of the return value 
# of hyperplane() (i.e. the vector d above) equals the 
# distance of both points.

X <- rbind(c(8,2),c(4,5))
hyperplane(X)
# If the order of the points is changed, the direction of d 
# (see details) may also change.



### ----<< Example 2 >>---- : unit vectors in R^3
X <- diag(rep(1,3))
hyperplane(X)
# The plane passing through all three unit vectors is given by 
#   -x - y - z + 1 = 0.
# These three points form a equilateral triangle on the plane 
# with side length sqrt(2) and hence area sqrt(3)/2.
# The norm of d (see details) equals twice this number. 
}

\keyword{multivariate}
\keyword{nonparametric}
