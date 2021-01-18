\name{Extract}
\alias{Extract.hyper2}
\alias{extract}
\alias{extractor}
\alias{[.hyper2}
\alias{[<-.hyper2}
\title{
Extract or replace parts of a hyper2 object
}
\description{
Extract or replace parts of a hyper2 object
}
\usage{
\method{[}{hyper2}(x, \dots)
\method{[}{hyper2}(x, index, \dots) <- value
}
\arguments{
  \item{x}{An object of class \code{hyper2}}
  \item{\dots}{Further arguments, currently ignored}
  \item{index}{A list with integer vector elements corresponding to the
      brackets whose power is to be replaced}
  \item{value}{Numeric vector of powers}
}
\details{
  These methods should work as expected, although the off-by-one issue
  might be a gotcha.

 For the extract method, a numeric vector is returned but the elements
 may be in any order (the brackets and indeed the powers are not stored
 in any particular order).

 The replace method, \code{H[L] <- value}, the index \code{L} is a list
 specifying the brackets to be overwritten; argument \code{value} is not
 recycled unless it is of length 1.

 If the index argument is missing, viz \code{H1[] <- H2}, this is a
 special case.  Argument \code{H1} must be a \code{hyper2} object, and
 the idiom effectively executes \code{H1[brackets(H2)] <- powers(H2)},
 but more efficiently (note that this operation is well-defined even
 though the order of the brackets is arbitrary).  This special case is
 included in the package because it has a very natural \code{C++}
 expression [function \code{overwrite()} in the \code{src/} directory]
 that was too neat to omit.

 Altering (incrementing or decrementing) the power of a single bracket
 is possible using idiom like \code{H[x] <- H[x] + 1}; this is
 documented at \code{Ops.hyper2}, specifically
 \code{hyper2_sum_numeric()}.
}

\value{
The extractor method returns a \code{hyper2} object, restricted to the
elements specified
}
\author{Robin K. S. Hankin}
\note{
Use \code{powers()} and \code{brackets()} to extract a numeric vector of
 powers or a list of integer vectors respectively.

Replacement idiom \code{H[x] <- val} cannot use non-trivial recycling.
This is because the elements of \code{H} are stored in an arbitrary
order, but the elements of \code{val} are stored in a particular order.
Also see function \code{hyper2_sum_numeric()}.

}
\seealso{\code{\link{hyper2}},\code{\link{Ops.hyper2}}}
\examples{
data(chess)

chess[1]
chess[c(1,2)]
chess[list(1,1:2)]
chess[1] <- 5

a <- hyper2(pnames=letters[1:8])
a[sapply(seq_len(8),seq_len)] <- 1
a

b <- a
b[1] <- 1000
b[2] <- 2000
b[] <- a      # b[1] overwritten; b[2] preserved


}
