\name{magnitude}
\alias{magnitude}
\alias{mod}
\alias{Mod}
\alias{Mod.clifford}
\title{Magnitude of a clifford object}
\usage{
\method{Mod}{clifford}(z)
}
\arguments{
  \item{z}{Clifford objects}
}
\description{

  Following Perwass, the magnitude of a multivector is defined as
  
  \deqn{\left|\left|A\right|\right| = \sqrt{A\ast A}}{omitted; see PDF}
  
  Where \eqn{A\ast A}{?} denotes the Euclidean scalar product
  \code{eucprod()}.  Recall that the Euclidean scalar product is never
  negative.
}
\author{Robin K. S. Hankin}
\note{

If you want the square, \eqn{\left|\left|A\right|\right|^2}{||A||^2} and
not \eqn{\left|\left|A\right|\right|}{||A||}, it is faster and more
accurate to use \code{eucprod(A)}, because this avoids a needless square
root.

There is a nice example of scalar product at \code{rcliff.Rd}.

}
\seealso{
  \code{\link{Ops.clifford}},
  \code{\link{Conj}},
  \code{\link{rcliff}}
}
\examples{


Mod(rcliff())


# Perwass, p68, asserts that if A is a k-blade, then (in his notation)
# AA == A*A.

# In package idiom, A*A == A \%star\% A:

A <- rcliff()          
Mod(A*A - A \%star\% A)  # meh

A <- rblade()
Mod(A*A - A \%star\% A)  # should be small



}

