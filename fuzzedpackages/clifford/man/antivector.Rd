\name{antivector}
\alias{antivector}
\alias{is.antivector}
\title{Antivectors or pseudovectors}
\description{
Antivectors or pseudovectors
}
\usage{
antivector(v, n = length(v))
is.antivector(C, include.pseudoscalar=FALSE)
}
\arguments{
  \item{v}{Numeric vector}
  \item{n}{Integer}
  \item{C}{Clifford object}
  \item{include.pseudoscalar}{Boolean: should the pseudoscalar be
    considered an antivector?}
}
\details{
An \dfn{antivector} is an \eqn{n}-dimensional Clifford object of all of
whose terms are of grade \eqn{n-1}.   

The pseudoscalar is a peculiar edge case.  Consider:

\preformatted{
  A <- clifford(list(c(1,2,3)))
  B <- A + clifford(list(c(1,2,4)))

> is.antivector(A)
[1] FALSE
> is.antivector(B)
[1] TRUE
> is.antivector(A,include.pseudoscalar=TRUE)
[1] TRUE
> is.antivector(B,include.pseudoscalar=TRUE)
[1] TRUE
}

One could argue that \code{A} should be an antivector as it is a term in
\code{B}, which is definitely an antivector.  Use
\code{include.pseudoscalar=TRUE} to ensure consistency in this case.
}

\note{
  An antivector is always a blade.
}
\references{

  Wikipedia contributors. (2018, July 20). \dQuote{Antivector}.  In
  \emph{Wikipedia, The Free Encyclopedia}. Retrieved 19:06, January 27,
  2020, from
  \url{https://en.wikipedia.org/w/index.php?title=Antivector&oldid=851094060}

}
\author{Robin K. S. Hankin}
\examples{
antivector(1:5)
}