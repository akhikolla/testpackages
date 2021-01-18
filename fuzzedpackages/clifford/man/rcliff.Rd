\name{rcliff}
\alias{rcliff}
\alias{rblade}
\title{Random clifford objects}
\description{Random Clifford algebra elements, intended as quick
  \dQuote{get you going} examples of \code{clifford} objects}
\usage{
rcliff(n=9, d=6, g=4, include.fewer=TRUE)
rblade(d=9, g=4)
}
\arguments{
  \item{n}{Number of terms}
  \item{d}{Dimensionality of underlying vector space}
  \item{g}{Maximum grade of any term}
  \item{include.fewer}{Boolean, with \code{FALSE} meaning to
    return a clifford object comprisig only terms of grade \code{g},
    and default \code{TRUE} meaning to include terms with
    grades less than \code{g}}
}
\details{

  Perwass gives the following lemma:

  Given blades \eqn{A_{\langle r\rangle}, B_{\langle
    s\rangle}, C_{\langle t\rangle}}{see PDF}, then

  \deqn{
    \langle
    A_{\langle r\rangle}
    B_{\langle s\rangle}
    C_{\langle t\rangle}
    \rangle_0
    =
    \langle
    C_{\langle t\rangle}
    A_{\langle r\rangle}
    B_{\langle s\rangle}
    \rangle_0
  }{
    omitted; see PDF
  }

  In the proof he notes in an intermediate step that 

    \deqn{
    \langle
    A_{\langle r\rangle}
    B_{\langle s\rangle}
    \rangle_t * 
    C_{\langle t\rangle}
    =
    C_{\langle t\rangle} * 
    \langle
    A_{\langle r\rangle}
    B_{\langle s\rangle}
    \rangle_t
    = 
    \langle
    C_{\langle t\rangle}
    A_{\langle r\rangle}
    B_{\langle s\rangle}
    \rangle_0.
  }{
    omitted; see PDF
  }

  Package idiom is shown in the examples.

  

}
\author{Robin K. S. Hankin}
\examples{

rcliff()
rcliff(d=3,g=2)
rcliff(3,10,7)
rcliff(3,10,7,include=TRUE)

x1 <- rcliff()
x2 <- rcliff()
x3 <- rcliff()

x1*(x2*x3) == (x1*x2)*x3  # should be TRUE


rblade()

# We can invert blades easily:
a <- rblade()
ainv <- rev(a)/scalprod(a)

zap(a*ainv)  # should be = 1
zap(ainv*a)  # should be = 1

# Perwass 2009, lemma 3.9:


A <- rblade(g=4)  # r=4
B <- rblade(g=5)  # s=5
C <- rblade(g=6)  # t=6

grade(A*B*C,0)-grade(C*A*B,0)   # geometric product uses '*'

# Intermediate step

x1 <- grade(A*B,7) \%star\% C
x2 <- C \%star\% grade(A*B,7)
x3 <- grade(C*A*B,0)

max(x1,x2,x3) - min(x1,x2,x3)   # should be small

}