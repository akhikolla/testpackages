\name{term}
\alias{term}
\alias{term}
\alias{terms}
\alias{is.term}
\alias{is.blade}
\alias{is.basisblade}
\title{Deal with terms}
\description{

By \dfn{basis vector}, I mean one of the basis vectors of the underlying
vector space \eqn{R^n}, that is, an element of the set
\eqn{\left\{e_1,\ldots,e_n\right\}}{{e_1,...,e_n}}.  A \dfn{term}
(sometimes a \dfn{basis blade} or \dfn{simple blade} is a wedge product
of basis vectors (or a geometric product of linearly independent basis
vectors), something like \eqn{e_{12}}{e_12} or \eqn{e_{12569}}{e_12569}.
  
From Perwass: a \dfn{blade} is the outer product of a number of
1-vectors (or, equivalently, the wedge product of linearly independent
1-vectors).  Thus \eqn{e_{12}=e_1\wedge e_2}{e_12=e_1 ^ e_2} and
\eqn{e_{12} + e_{13}=e_1\wedge(e_2+e_3)}{e_12+e_13=e1^(e2+e3)} are
blades, but \eqn{e_{12} + e_{34}}{e_12+e_34} is not.

Function \code{is.blade()} is not currently implemented: there is no
easy way to detect whether a Clifford object is a product of 1-vectors.

}

\usage{
terms(x)
is.blade(x)
is.basisblade(x)
}
\arguments{
  \item{x}{Object of class \code{clifford}}
}
\details{

  \itemize{
    
    \item Functions \code{terms()} and \code{coeffs()} are the
    extraction methods.  These are unordered vectors but the ordering is
    consistent between them (an extended discussion of this phenomenon
    is presented in the \code{mvp} package).

    \item Function \code{term()} returns clifford object that comprises
    a single term with unit coefficient.
    
    \item Function \code{is.basisterm()} returns \code{TRUE} if its
    argument has only a single term, or is a nonzero scalar; the zero
    clifford object is not considered to be a basis term.

  }
}
\references{
  C. Perwass. \dQuote{Geometric algebra with applications in
  engineering}.   Springer, 2009.
}
\author{Robin K. S. Hankin}
\seealso{\code{\link{clifford}}}
\examples{

x <- rcliff()
terms(x)

is.basisblade(x)


a <- as.1vector(1:3)
b <- as.1vector(c(0,0,0,12,13))

a \%^\% b # a blade

}
