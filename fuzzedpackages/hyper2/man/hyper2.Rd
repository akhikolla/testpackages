\name{hyper2}
\alias{hyper2}
\alias{brackets}
\alias{brackets.hyper2}
\alias{powers}
\alias{powers.hyper2}
\alias{pnames}
\alias{pnames.hyper2}
\alias{pnames<-}
\alias{pnames<-.hyper2}
\alias{size}
\alias{size.hyper2}
\alias{as.hyper2}
\alias{is_valid_hyper2}
\alias{is_constant}
\title{
Basic functions in the hyper2 package
}
\description{
Basic functions in the hyper2 package
}
\usage{
hyper2(L, d, pnames = NA)
\method{brackets}{hyper2}(H)
\method{powers}{hyper2}(H)
\method{pnames}{hyper2}(H)
size(H)
as.hyper2(L,d,pnames)
is_valid_hyper2(L,d,pnames)
is_constant(H)
}
\arguments{
  \item{H}{A hyper2 object}
  \item{L}{A list whose elements specify the brackets of a hyper2 object}
  \item{d}{A vector of powers; \code{hyper2()} recycles \emph{only if}
    \code{d} is of length 1}
  \item{pnames}{
    A character vector specifying the names of \eqn{p_1} through \eqn{p_n}.}
}
\details{
These are the basic functions of the hyper2 package.  Function
\code{hyper()} is the low-level creator function; \code{as.hyper2()} is
a bit more user-friendly and attempts to coerce its arguments into a
suitable form; for example, a matrix is interpreted as rows of brackets.

Function \code{is_valid_hyper2()} tests for valid input, returning a
Boolean.

Function \code{size()} returns the (nominal) length \eqn{n} of
nonnegative vector \eqn{p=\left(p_1,\ldots,p_n\right)}{p=(p_1,...,p_n)}
where \eqn{p_1+\cdots+p_n=1}{p_1+...+p_n=1}.  Note that function
\code{lhyper2()} takes the vector
\eqn{p=\left(p_1,\ldots,p_{n-1}\right)}{p=(p_1,...,p_{n-1})}.

}
\author{Robin K. S. Hankin}
\seealso{\code{\link{Ops.hyper2}},\code{\link{Extract}},\code{\link{loglik}}}
\examples{
hyper2(list(1,2,3,1:2,2:3),1:5)
}
