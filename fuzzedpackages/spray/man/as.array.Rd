\name{as.array}
\alias{as.array}
\alias{as.array.spray}
\alias{dim.spray}
\title{Coerce spray objects to arrays
}
\description{
Coerces spray objects to arrays.  Includes off-by-one functionality via
option \code{offbyone}.
}
\usage{
\method{as.array}{spray}(x, offbyone=FALSE, compact=FALSE, ...)
\method{dim}{spray}(x)
}
\arguments{
  \item{x}{spray object}
  \item{offbyone}{Boolean with default \code{FALSE} meaning to interpret
    the index entries as positions in their dimension, and \code{TRUE}
    meaning to add one to index values so that zero entries appear in
    the first place}
  \item{compact}{Boolean with default \code{FALSE} meaning to translate
    the spray as is, and \code{TRUE} meaning to add constants to each
    column of the index matrix so that the resulting array is as small
    as possible}
  \item{...}{Further arguments, currently ignored}
}
\details{
  
  Argument \code{offbyone} defaults to \code{FALSE}; but if it is set to
  \code{TRUE}, it effectively adds one from the index matrix, so a zero
  entry in the index matrix means the first position in that dimension.

  After the subtraction, if performed, the function will not operate if
  any index is less than 1.
  
}
\value{

  Returns an array of dimension \code{dim(S)}.  The \dQuote{meat} of the
  function is

\preformatted{
    out <- array(0, dS)
    out[ind] <- value(S)
}
}
\author{
Robin K. S. Hankin
}

\examples{

M <- matrix(sample(0:4,28,replace=TRUE),ncol=4)
S <- spray(M,sample(7),addrepeats=TRUE)
A <- as.array(S,offbyone=TRUE) 



S <- spray(matrix(sample(1:4,28,replace=TRUE),ncol=4),sample(7))
A <- as.array(S)   # S has no zero indices

stopifnot(all(S[index(S),drop=TRUE] == A[index(S)])) 
}
\keyword{symbolmath}

