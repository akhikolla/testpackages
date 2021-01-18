\name{Ops.clifford}
\alias{\%.\%}
\alias{\%^\%}
\alias{\%X\%}
\alias{\%star\%}
\alias{\% \%}
\alias{\%euc\%}
\alias{\%o\%}
\alias{\%_|\%}
\alias{\%|_\%}
\alias{\%.\%.clifford}
\alias{\%^\%.clifford}
\alias{\%X\%.clifford}
\alias{\%star\%.clifford}
\alias{\% \%.clifford}
\alias{\%euc\%.clifford}
\alias{\%o\%.clifford}
\alias{\%_|\%.clifford}
\alias{\%|_\%.clifford}
\alias{Ops.clifford}
\alias{Ops}
\alias{clifford_negative}
\alias{clifford_times_clifford}
\alias{geoprod}
\alias{geometric_prod}
\alias{geometric_product}
\alias{clifford_times_scalar}
\alias{clifford_plus_clifford}
\alias{clifford_plus_numeric}
\alias{clifford_plus_scalar}
\alias{clifford_star_clifford}
\alias{star}
\alias{clifford_cross_clifford}
\alias{cross}
\alias{clifford_fatdot_clifford}
\alias{fatdot}
\alias{clifford_lefttick_clifford}
\alias{lefttick}
\alias{clifford_righttick_clifford}
\alias{righttick}
\alias{clifford_power_scalar}
\alias{clifford_eq_clifford}
\alias{clifford_inverse}
\alias{dot}
\alias{wedge}
\alias{clifford_dot_clifford}
\alias{cliffdotprod}
\alias{clifford_wedge_clifford}
\alias{wedge}
\alias{maxyterm}
\alias{scalprod}
\alias{scalar_product}
\alias{eucprod}
\alias{euclid\_product}
\alias{euclidean\_product}
\title{Arithmetic Ops Group Methods for \code{clifford} objects}
\description{
  Allows arithmetic operators to be used for
  multivariate polynomials  such as addition, multiplication,
  integer powers, etc.
}
\usage{
\method{Ops}{clifford}(e1, e2)
clifford_negative(C)
geoprod(C1,C2)
clifford_times_scalar(C,x)
clifford_plus_clifford(C1,C2)
clifford_eq_clifford(C1,C2)
clifford_inverse(C)
cliffdotprod(C1,C2)
fatdot(C1,C2)
lefttick(C1,C2)
righttick(C1,C2)
wedge(C1,C2)
scalprod(C1,C2=rev(C1),drop=TRUE)
eucprod(C1,C2=C1,drop=TRUE)
maxyterm(C1,C2=as.clifford(0))
C1 \%.\% C2
C1 \%^\% C2
C1 \%X\% C2
C1 \%star\% C2
C1 \% \% C2
C1 \%euc\% C2
C1 \%o\% C2
C1 \%_|\% C2
C1 \%|_\% C2
}
\arguments{
  \item{e1,e2,C,C1,C2}{Objects of class \code{clifford}}
  \item{x}{Scalar, length one numeric vector}
  \item{drop}{Boolean, with default \code{TRUE} meaning to return the
    constant coerced to numeric, and \code{FALSE} meaning to return a
    (constant) Clifford object}
}
\details{

  The function \code{Ops.clifford()} passes unary and binary arithmetic
  operators \dQuote{\code{+}}, \dQuote{\code{-}}, \dQuote{\code{*}},
  \dQuote{\code{/}} and \dQuote{\code{^}} to the appropriate specialist
  function.

  Functions \code{c_foo()} are low-level helper functions that wrap the
  \code{C} code; function \code{maxyterm()} returns the maximum index
  in the terms of its arguments.

  The package has several binary operators:

  \tabular{lll}{\cr
    Geometric product \tab \code{A*B = geoprod(A,B)} \tab
    \eqn{\displaystyle AB=\sum_{r,s}\left\langle A\right\rangle_r\left\langle
      B\right\rangle_s}{see PDF}\cr
    Inner product\tab \code{A \%.\% B = cliffdotprod(A,B)}\tab
    \eqn{\displaystyle A\cdot B=\sum_{r\neq 0\atop s\ne 0}\left\langle\left\langle A\right\rangle_r\left\langle
      B\right\rangle_s\right\rangle_{\left|s-r\right|}}{see PDF}\cr
    Outer product\tab\code{A \%^\% B = wedge(A,B)}\tab
    \eqn{\displaystyle A\wedge B=\sum_{r,s}\left\langle\left\langle A\right\rangle_r\left\langle
      B\right\rangle_s\right\rangle_{s+r}}{see PDF}\cr
    Fat dot product\tab\code{A \%o\% B = fatdot(A,B)}\tab
    \eqn{\displaystyle A\bullet B=\sum_{r,s}\left\langle\left\langle A\right\rangle_r\left\langle
      B\right\rangle_s\right\rangle_{\left|s-r\right|}}{see PDF}\cr
    Left contraction\tab\code{A \%_|\% B = lefttick(A,B)}\tab
    \eqn{\displaystyle A\rfloor B=\sum_{r,s}\left\langle\left\langle A\right\rangle_r\left\langle
      B\right\rangle_s\right\rangle_{s-r}}{see PDF}\cr
    Right contraction\tab\code{A \%|_\% B = righttick(A,B)}\tab
    \eqn{\displaystyle A\lfloor B=\sum_{r,s}\left\langle\left\langle A\right\rangle_r\left\langle
      B\right\rangle_s\right\rangle_{r-s}}{see PDF}\cr
    Cross product\tab\code{A \%X\% B = cross(A,B)}\tab
    \eqn{\displaystyle A\times
    B=\frac{1}{2_{\vphantom{j}}}\left(AB-BA\right)}{see PDF}\cr
    Scalar product\tab\code{A \%star\% B = star(A,B)}\tab
    \eqn{\displaystyle A\ast B=\sum_{r,s}\left\langle\left\langle A\right\rangle_r\left\langle
      B\right\rangle_s\right\rangle_0}{see PDF}\cr
    Euclidean product\tab\code{A \%euc\% B = eucprod(A,B)}\tab
    \eqn{\displaystyle A\star B= A\ast B^\dag}{see PDF}\cr
  }

  In R idiom, the geometric product \code{geoprod(.,.)} has to be
  indicated with a \dQuote{\code{*}} (as in \code{A*B}) and so the
  binary operator must be \code{\%*\%}: we need a different idiom for
  scalar product, which is why \code{\%star\%} is used.

  Because geometric product is often denoted by juxtaposition, package
  idiom includes \code{a \% \% b} for geometric product.

  Function \code{clifford_inverse()} is problematic as nonnull blades
  always have an inverse; but function \code{is.blade()} is not yet
  implemented.  Blades (including null blades) have a pseudoinverse, but
  this is not implemented yet either.
  
  The \dfn{scalar product} of two clifford objects is defined as the
  zero-grade component of their geometric product:
  
  \deqn{
    A\ast B=\left<AB\right>_0\qquad{\mbox{NB: notation used by both Perwass
    and Hestenes}}
  }{
    omitted; see PDF
  }

  In package idiom the scalar product is given by \code{A \%star\% B} or
  \code{scalprod(A,B)}.  Hestenes and Perwass both use an asterisk for
  scalar product as in \dQuote{\eqn{A*B}}, but in package idiom, the
  asterisk is reserved for geometric product.

  \strong{Note: in the package, \code{A*B} is the geometric product}.
  
  The \dfn{Euclidean product} (or \dfn{Euclidean scalar product}) of two
  clifford objects is defined as 

  \deqn{
    A\star B=A\ast B^\dag=\left<AB^\dag\right>_0\qquad{\mbox{Perwass}}
  }{
    omitted: see PDF
  }

  where \eqn{B^\dag}{?} denotes Conjugate [as in \code{Conj(a)}].  In
  package idiom the Euclidean scalar product is given by
  \code{eucprod(A,B)} or \code{A \%euc\% B}, both of which return
  \code{A * Conj(B)}.

  Note that the scalar product \eqn{A\ast A} can be positive or negative
  [that is, \code{A \%star\% A} may be any sign], but the Euclidean
  product is guaranteed to be non-negative [that is, \code{A \%euc\% A}
  is always positive or zero].

  Dorst defines the left and right contraction (Chisholm calls these the
  left and right inner product) as \eqn{A\rfloor B}{see PDF} and
  \eqn{A\lfloor B}.  See the vignette for more details.

}
\value{

  The high-level functions documented here return an object of
  \code{clifford}.  But don't use the low-level functions.

}
\author{Robin K. S. Hankin}
\seealso{\code{\link{scalprod}}}
\examples{

u <- rcliff(5)
v <- rcliff(5)
w <- rcliff(5)

u*v

u^3


u+(v+w) == (u+v)+w              # should be TRUE 
u*(v*w) == (u*v)*w              # should be TRUE 
u \%^\% v == (u*v-v*u)/2        # should be TRUE 

# Now if x,y,z are _vectors_ we would have:

x <- as.1vector(5)
y <- as.1vector(5)
x*y == x\%.\%y + x\%^\%y        # should be TRUE
x \%^\% y == x \%^\% (y + 3*x)  # should be TRUE 

## Inner product is "%.%" is not associative:
 rcliff(5,g=2) -> x
 rcliff(5,g=2) -> y
 rcliff(5,g=2) -> z
x \%.\% (y \%.\% z)
(x \%.\% y) \%.\% z

## Geometric product *is* associative:
x * (y * z)
(x * y) * z


}

