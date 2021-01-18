\name{Extract.clifford}
\alias{[.clifford}
\alias{extract.clifford}
\alias{extract}
\alias{replace.clifford}
\alias{replace}
\alias{[<-.clifford}
\alias{coeffs}
\alias{coeffs<-}
\alias{coeffs<-.clifford}
\title{Extract or Replace Parts of a clifford}
\description{
  Extract or replace subsets of cliffords.
}
\usage{
\method{[}{clifford}(C, index, \dots)
\method{[}{clifford}(C, index, \dots) <- value
coeffs(x)
coeffs(x) <- value
}
\arguments{
  \item{C,x}{A clifford object}
  \item{index}{elements to extract or replace}
  \item{value}{replacement value}
  \item{...}{Further arguments}
}
\details{

Extraction and replacement methods.  The extraction method uses
\code{getcoeffs()} and the replacement method uses low-level helper
function \code{c_overwrite()}.

In the extraction function \code{a[index]}, if \code{index} is a list,
further arguments are ignored.  If not, the dots are used.

Replacement methods using list-valued \code{index}, as in \code{A[i] <-
value} uses an ugly hack if \code{value} is zero.

Idiom such as \code{a[] <- b} follows the \code{spray} package.  If
\code{b} is a length-one scalar, then \code{coeffs(a) <- b} has the same
effect as \code{a[] <- b}.

Functions \code{terms()} [see \code{term.Rd}]and \code{coeffs()} are the
extraction methods.  These are unordered vectors but the ordering is
consistent between them (an extended discussion of this phenomenon is
presented in the \code{mvp} package).

Function \code{coeffs<-()} (idiom \code{coeffs(a) <- b}) sets all
coefficients of \code{a} to \code{b}.  This has the same effect as
\code{a[] <- b}.
 
}
\seealso{\code{\link{Ops.clifford}},\code{\link{clifford}},\code{\link{term}}}
\examples{
A <- clifford(list(1,1:2,1:3),1:3)
B <- clifford(list(1:2,1:6),c(44,45))

A[1,c(1,3,4)]

A[2:3, 4] <- 99
A[] <- B

}


