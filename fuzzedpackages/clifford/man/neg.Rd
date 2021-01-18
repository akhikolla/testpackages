\name{neg}
\alias{neg}
\title{Grade negation}
\description{
The grade \eqn{r} negation operation applied to Clifford multivector
\eqn{A} changes the sign of the grade \eqn{r} component of \eqn{A}.  It
is formally defined as \eqn{A-2\left<A\right>_r}{A-2<A>_r}.
}
\usage{
neg(C, n)
}
\arguments{
  \item{C}{Clifford object}
  \item{n}{Integer vector indicating grades to negate}
}
\details{
The function is algebraically equivalent to
\code{function(C,n){C-2*grade(C,n)}} but uses faster and more efficient idiom.
}
\references{
A. Acus and A. Dargys 2018.  \dQuote{The inverse of a multivector:
  beyond the threshold \eqn{p+q=5}}.  \emph{Advances in applied Clifford
  Algebras}, 28:65
}
\author{Robin K. S. Hankin}
\seealso{\code{\link{Conj}}}
\examples{

A <- rcliff()
neg(A,1:2) == A-grade(A,1:2)  # should be TRUE

}