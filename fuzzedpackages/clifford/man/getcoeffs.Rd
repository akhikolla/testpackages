\name{getcoeffs}
\alias{getcoeffs}
\alias{const}
\alias{constant}
\alias{const<-}
\alias{const<-.clifford}
\alias{constant<-}
\alias{constant<-.clifford}
\title{Get coefficients of a Clifford object}
\description{
Access specific coefficents of a Clifford object using a list of terms.
}
\usage{
getcoeffs(C, B)
const(C,drop=TRUE)
\method{const}{clifford}(x) <- value
}
\arguments{
  \item{C,x}{Clifford object}
  \item{B}{List of terms}
  \item{value}{Replacement value}
  \item{drop}{Boolean, with default \code{TRUE} meaning to return the
    constant coerced to numeric, and \code{FALSE} meaning to return a
    (constant) Clifford object}
}
\details{

Extractor method for specific terms.  Function \code{const()} returns
the constant element of a Clifford object.  Note that \code{const(C)}
returns the same as \code{grade(C,0)}, but is faster.

The slightly awkward R idiom in \code{const<-()} is to ensure numerical
accuracy; see examples.

}
\author{Robin K. S. Hankin}
\seealso{\code{\link{clifford}}}
\examples{

X <- clifford(list(1,1:2,1:3,3:5),6:9)

getcoeffs(X,1:2)

X <- X+1e300

const(X) # should be 1e300

const(X) <- 0.6
const(X) # should be 0.6, no numerical error

}
