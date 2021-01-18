\name{spray_cpp}
\alias{spray_rcpp}
\alias{spray_accessor}
\alias{spray_add}
\alias{spray_crush}
\alias{spray_maker}
\alias{spray_mult}
\alias{spray_overwrite}
\alias{spray_setter}
\alias{spray_equality}
\alias{spray_asum_include}
\alias{spray_asum_exclude}
\alias{spray_deriv}
\alias{spray_pmax}
\alias{spray_pmin}
\alias{spray_power}
\alias{spray_spray_accessor}
\alias{spray_spray_add}
\alias{spray_spray_asum_exclude}
\alias{spray_spray_asum_include}
\alias{spray_spray_crush}
\alias{spray_spray_deriv}
\alias{spray_spray_equality}
\alias{spray_spray_maker}
\alias{spray_spray_mult}
\alias{spray_spray_overwrite}
\alias{spray_spray_pmax}
\alias{spray_spray_pmin}
\alias{spray_spray_setter}
\alias{spray_spray_power}
\title{
Low-level functions that call C++ source code
}
\description{
Low-level functions that call C++ source code, as detailed in the
automatically generated \code{RcppExports.R} file.
}
\usage{
spray_maker(M, d)
spray_add(M1, d1, M2, d2)
spray_mult(M1, d1, M2, d2)
spray_overwrite(M1, d1, M2, d2)
spray_accessor(M, d, Mindex)
spray_setter(M1, d1, M2, d2)
spray_equality(M1, d1, M2, d2)
spray_asum_include(M,d,n)
spray_asum_exclude(M,d,n)
spray_deriv(M,d,n)
spray_pmax(M1,d1,M2,d2)
spray_pmin(M1,d1,M2,d2)
spray_power(M,d,pow)
spray_spray_accessor()
spray_spray_add()
spray_spray_asum_exclude()
spray_spray_asum_include()
spray_spray_deriv()
spray_spray_equality()
spray_spray_maker()
spray_spray_mult()
spray_spray_overwrite()
spray_spray_pmax()
spray_spray_pmin()
spray_spray_setter()
spray_spray_power()
}
\arguments{
  \item{M,M1,M2,Mindex}{Integer valued matrices with rows corresponding to
    array indices}
  \item{d,d1,d2}{Vector of values corresponding to nonzero array
    entries}
  \item{n}{Integer vector corresponding to dimensions to sum over for
    the sum functions}
  \item{pow}{Nonnegative integer for \code{spray_power()}}
}
\value{
These functions return a two-element list which is coerced to
an object of class \code{spray}  by function \code{spraymaker()}.
}
\author{Robin K. S. Hankin}
\note{
  These functions aren't really designed for the end-user.

  Function \code{spray_equality()} cannot simply check for equality of
  \code{$value} because the order of the index rows is not specified in
  a spray object.  Function \code{spray_crush()} has been removed as it
  is redundant.
}
\seealso{
\code{\link{spraymaker}},\code{\link{spray}}
}
\keyword{symbolmath}
