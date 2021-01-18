\name{drop}
\alias{drop}
\title{Drop redundant information}
\description{Coerce constant Clifford objects to numeric}
\usage{drop(C)}
\arguments{\item{C}{Clifford object}}
\details{If its argument is a constant clifford object, coerce to numeric.}
\author{Robin K. S. Hankin}
\seealso{\code{\link{grade}},\code{\link{getcoeffs}}}
\examples{
drop(as.clifford(5))

const(rcliff())
const(rcliff(),drop=FALSE)

}
