\name{rfalg}
\alias{rfalg}
\alias{rfreealg}
\title{Random free algebra objects}
\description{Random elements of the free algebra, intended as quick
  \dQuote{get you going} examples of \code{freealg} objects}
\usage{
rfalg(n=7, distinct=3, maxsize=4, include.negative=FALSE)
}
\arguments{
  \item{n}{Number of terms to generate}
  \item{distinct}{Number of distinct symbols to use}
  \item{maxsize}{Maximum number of symbols in any word}
  \item{include.negative}{Boolean, with default \code{FALSE} meaning to
    use only positive symbols (lower-case letters) and \code{TRUE}
    meaning to use upper-case letters as well, corresponding to the
    inverse of the lower-case symbols}
}
\details{
What you see is what you get, basically.  A term such as
\code{aabaAbBB} will be simplified to \code{aabbBB}.
}
\author{Robin K. S. Hankin}
\examples{
rfalg()
rfalg()^3

constant(rfalg())

}