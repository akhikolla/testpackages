\name{zap}
\alias{zap}
\alias{zapsmall}
\alias{zaptiny}
\alias{zap.clifford}
\title{Zap small values in a clifford object}
\description{Generic version of \code{zapsmall()}}
\usage{
zap(x, drop=TRUE, digits = getOption("digits"))
}
\arguments{
  \item{x}{Clifford object}
  \item{drop}{Boolean with default \code{TRUE} meaning to 
    coerce the output to numeric with \code{drop()}}
  \item{digits}{number of digits to retain}
}
\details{
  
  Given a clifford object, coefficients close to zero are
  \sQuote{zapped}, i.e., replaced by \sQuote{0} in much the same way as
  \code{base::zapsmall()}.

  The function should be called \code{zapsmall()}, and dispatch to the
  appropriate base function, but I could not figure out how to do this
  with S3 (the docs were singularly unhelpful) and gave up. 

  Note, this function actually changes the numeric value, it is not just
  a print method.
}
\author{Robin K. S. Hankin}
\examples{

a <- clifford(sapply(1:10,seq_len),90^-(1:10))
zap(a)
options(digits=3)
zap(a)


a-zap(a)   # nonzero

B <- rblade(g=3)
mB <- B*rev(B)
zap(mB)
drop(mB)
}
