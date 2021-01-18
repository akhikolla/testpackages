\name{lowlevel}
\alias{c_identity}
\alias{c_add}
\alias{c_multiply}
\alias{c_power}
\alias{c_equal}
\alias{c_grade}
\alias{c_getcoeffs}
\alias{c_innerprod}
\alias{c_outerprod}
\alias{c_fatdotprod}
\alias{c_lefttickprod}
\alias{c_righttickprod}
\alias{c_overwrite}
\title{Low-level helper functions for \code{clifford} objects}
\description{
  Helper functions for \code{clifford} objects, written in \code{C} using the
  \code{STL} map class.
}
\usage{
c_identity(L, p, m)
c_grade(L, c, m, n)
c_add(L1, c1, L2, c2, m)
c_multiply(L1, c1, L2, c2, m, sig)
c_power(L, c, m, p, sig)
c_equal(L1, c1, L2, c2, m)
c_overwrite(L1, c1, L2, c2, m)
}
\arguments{
  \item{L,L1,L2}{Lists of terms}
  \item{c1,c2,c}{Numeric vectors of coefficients}
  \item{m}{Maximum entry of terms}
  \item{n}{Grade to extract}
  \item{p}{Integer power}
  \item{sig}{Positive integer representing number of \eqn{+1} on main
    diagonal of quadratic form}
}
\details{

  The functions documented here are low-level helper functions that wrap
  the \code{C} code.  They are called by functions like
  \code{clifford_plus_clifford()}, which are themselves called by the
  binary operators documented at \code{Ops.clifford.Rd}.

  Function \code{clifford_inverse()} is problematic as nonnull blades
  always have an inverse; but function \code{is.blade()} is not yet
  implemented.  Blades (including null blades) have a pseudoinverse, but
  this is not implemented yet either.

}
\value{

  The high-level functions documented here return an object of
  \code{clifford}.  But don't use the low-level functions.

}
\author{Robin K. S. Hankin}
\seealso{\code{\link{Ops.clifford}}}