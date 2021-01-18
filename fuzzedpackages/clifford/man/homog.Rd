\name{homog}
\alias{homog}
\alias{is.homog}
\alias{is.homogenous}
\alias{homogenous}
\title{Homogenous Clifford objects}
\description{

A clifford object is homogenous if all its terms are the same grade.  A
scalar (including the zero clifford object) is considered to be
homogenous.   This ensures that

\code{is.homog(grade(C,n))} always returns \code{TRUE}.

}
\usage{
is.homog(C)
}
\arguments{\item{C}{Object of class clifford}}
\details{
  Homogenous clifford objects have a multiplicative inverse.
}
\author{Robin K. S. Hankin}
\examples{
is.homog(rcliff())
is.homog(rcliff(include.fewer=FALSE))
}