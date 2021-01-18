\name{even}
\alias{even}
\alias{is.even}
\alias{is.odd}
\alias{evenpart}
\alias{oddpart}
\title{Even and odd clifford objects}
\description{
A clifford object is \emph{even} if every term has even grade, and
\emph{odd} if every term has odd grade.

Functions \code{is.even()} and \code{is.odd()} test a clifford object
for evenness or oddness.

Functions \code{evenpart()} and \code{oddpart()} extract the even 
or odd terms from a clifford object, and we write \eqn{A_+}{A+} and
\eqn{A_-}{A-} respectively; we have \eqn{A=A_+ + A_-}

}
\usage{
is.even(C)
is.odd(C)
evenpart(C)
oddpart(C)
}
\arguments{\item{C}{Clifford object}}
\author{Robin K. S. Hankin}
\seealso{\code{\link{grade}}}
\examples{

A <- rcliff()
A == evenpart(A) + oddpart(A) # should be true

}