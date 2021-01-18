\name{Ops.freealg}
\alias{Ops.freealg}
\alias{ops}
\alias{free_negative}
\alias{freealg_negative}
\alias{free_times_free}
\alias{free_times_scalar}
\alias{free_plus_free}
\alias{free_plus_numeric}
\alias{free_power_scalar}
\alias{free_eq_free}
\alias{free_equal_free}
\alias{lowlevel_simplify}
\alias{lowlevel_free_prod}
\alias{lowlevel_free_sum}
\alias{lowlevel_free_power}
\title{Arithmetic Ops methods for the the free algebra}
\description{
  Arithmetic operators for manipulation of freealg objects
  such as addition, multiplication, powers, etc
}
\usage{
\method{Ops}{freealg}(e1, e2)
free_negative(S)
free_power_scalar(S,n)
free_eq_free(e1,e2)
free_plus_numeric(S,x)
free_plus_free(e1,e2)
lowlevel_simplify(words,coeffs)
lowlevel_free_prod(words1,coeffs1,words2,coeffs2)
lowlevel_free_sum(words1,coeffs1,words2,coeffs2)
lowlevel_free_power(words,coeffs,n)
}
\arguments{
  \item{S,e1,e2}{Objects of class \code{freealg}}
  \item{n}{An integer, possibly non-positive}
  \item{x}{Scalar value}
  \item{words,words1,words2}{A list of words, that is, a list of integer vectors
    representing the variables in each term}
  \item{coeffs,coeffs1,coeffs2}{Numeric vector representing the
    coefficients of each word}
}
\details{

  The function \code{Ops.freealg()} passes binary arithmetic operators
  (\dQuote{\code{+}}, \dQuote{\code{-}}, \dQuote{\code{*}},
  \dQuote{\code{^}}, and \dQuote{\code{==}}) to the appropriate
  specialist function.

  The caret, as in \code{a^n}, denotes arithmetic exponentiation, as in
  \code{x^3==x*x*x}.

  Functions \code{lowlevel_foo()} are low-level functions that interface
  directly with the \code{C} routines in the \code{src/} directory and
  are not intended for the end-user.
}
\author{Robin K. S. Hankin}
\examples{
rfalg()
as.freealg("1+x+xy+yx")  # variables are non-commutative
as.freealg("x") * as.freealg("X") # upper-case letters are lower-case inverses



constant(as.freealg("x+y+X+Y")^6)  # OEIS sequence A035610
}
