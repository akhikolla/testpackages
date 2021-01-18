\name{numeric_to_clifford}
\alias{numeric_to_clifford}
\alias{as.cliffvector}
\alias{as.1vector}
\alias{is.1vector}
\alias{is.scalar}
\alias{as.scalar}
\alias{scalar}
\alias{basis}
\alias{e}
\alias{pseudoscalar}
\alias{as.pseudoscalar}
\alias{is.pseudoscalar}
\title{Coercion from numeric to Clifford form}
\description{
Given a numeric value or vector, return a Clifford algebra element}
\usage{
numeric_to_clifford(x)
as.1vector(x)
is.1vector(x)
scalar(x=1)
as.scalar(x=1)
is.scalar(C)
basis(n,x=1)
e(n,x=1)
pseudoscalar(n,x=1)
as.pseudoscalar(n,x=1)
is.pseudoscalar(C)
}
\arguments{
  \item{x}{Numeric vector}
  \item{n}{Integer specifying dimensionality of underlying vector space} 
  \item{C}{Object possibly of class Clifford}
}
\details{

  Function \code{as.scalar()} takes a length-one numeric vector and
  returns a Clifford scalar of that value (to extract the scalar
  component of a multivector, use \code{const()}).

  Function \code{as.1vector()} takes a numeric vector and returns the
  linear sum of length-one blades with coefficients given by \code{x};
  function \code{is.1vector()} returns \code{TRUE} if every term is of
  grade 1.

  Function \code{pseudoscalar(n)} returns a pseudoscalar of
  dimensionality \code{n} and function \code{is.pseudoscalar()} checks
  for a Clifford object being a pseudoscalar.
  
  Function \code{numeric_to_vector()} dispatches to either
  \code{as.scalar()} for length-one vectors or \code{as.1vector()}
  if the length is greater than one.

  Function \code{basis()} returns a clifford element comprising of
  single blade with grade 1; function \code{e()} is a synonym.

}
\author{Robin K. S. Hankin}
\seealso{\code{\link{getcoeffs}}}
\examples{

as.scalar(6)
as.1vector(1:8)

Reduce(`+`,sapply(seq_len(7),function(n){e(seq_len(n))},simplify=FALSE))

pseudoscalar(6)

pseudoscalar(7,5) == 5*pseudoscalar(7)  # should be true


}

