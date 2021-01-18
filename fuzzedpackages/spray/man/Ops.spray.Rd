\name{Ops.spray}
\alias{Ops.spray}
\alias{Ops}
\alias{spray_negative}
\alias{spray_times_spray}
\alias{spray_times_scalar}
\alias{spray_plus_spray}
\alias{spray_plus_scalar}
\alias{spray_power_scalar}
\alias{spray_eq_spray}
\title{Arithmetic Ops Group Methods for sprays}
\description{
  Allows arithmetic operators to be used for spray calculations, 
  such as addition, multiplication, division, integer powers, etc.
  Objects of class spray are interpreted as sparse multivariate polynomials.
}
\usage{
\method{Ops}{spray}(e1, e2 = NULL)
spray_negative(S)
spray_times_spray(S1,S2)
spray_times_scalar(S,x)
spray_plus_spray(S1,S2)
spray_plus_scalar(S,x)
spray_power_scalar(S,n)
spray_eq_spray(S1,S2)
}
\arguments{
  \item{e1,e2,S,S1,S2}{Objects of class spray, here interpreted as
    sparse multivariate polynomials}
  \item{x}{Real valued scalar}
  \item{n}{Non-negative integer}
}
\details{

  The function \code{Ops.spray()} passes unary and binary arithmetic
  operators (\dQuote{\code{+}}, \dQuote{\code{-}}, \dQuote{\code{*}},
  \dQuote{\code{/}},\dQuote{\code{==}}, and \dQuote{\code{^}}) to the
  appropriate specialist function.

  The most interesting operators are \dQuote{\code{*}} and
  \dQuote{\code{+}} which execute multivariate polynomial multiplication
  and addition respectively.

  Testing for equality uses \code{spray_eq_spray()}.  Note that
  \code{spray_eq_spray(S1,S2)} is algebraically equivalent to
  \code{is.zero(S1-S2)}, but faster (\code{FALSE} is returned as soon as
  a mismatch is found).

}
\value{
The functions all return spray objects except  \dQuote{\code{==}}, which
returns a logical.
}
\author{
Robin K. S. Hankin
}
\note{
Notes here
}
\seealso{
\code{\link{ooom}}
}
\examples{


M <- matrix(sample(0:3,21,replace=TRUE),ncol=3)
a <- spray(M,sample(7))
b <- homog(3,4)


# arithmetic operators mostly work as expected:
a + 2*b
a - a*b^2/4
a+b

S1 <- spray(partitions::compositions(4,3))
S2 <- spray(diag(3))  # S2 = x+y+z


stopifnot( (S1+S2)^3 == S1^3 + 3*S1^2*S2 + 3*S1*S2^2 + S2^3 )






}
\keyword{symbolmath}
