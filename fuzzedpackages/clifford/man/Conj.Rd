\name{Conj}
\alias{Conj.clifford}
\alias{Conj}
\alias{dagger}
\title{Conjugate of a Clifford object}
\description{The \dQuote{conjugate} of a Clifford object is defined by
  Perwass in definition 2.9, p59.}
\usage{\method{Conj}{clifford}(z)}
\arguments{\item{z}{Clifford object}}
\details{

  Perwass uses a dagger to indicate Conjugates, as in \eqn{A^\dag}.  If

  \deqn{
    A_{\langle k\rangle}
    =
    \bigwedge_{i=1}^k a_i
  }{
    omitted
  }

  Then
    
  \deqn{
    A_{\langle k\rangle}^\dag
    =
    \left(a_1\wedge\ldots\wedge a_k\right)^\dag
    =
    a_k^\dag\wedge\ldots\wedge a_1^\dag
    =
    \bigwedge_{i=1}^k a_i^\dag
  }{
    omitted
  }

  He  gives the following theorem (3.58, p70]:
  
  Given blades \eqn{A_{\langle k\rangle}, B_{\langle l\rangle}}{see
    PDF}, then
  
  \deqn{
    \left(
    A_{\langle k\rangle}\wedge B_{\langle l\rangle}
    \right)^\dag
    = 
    B_{\langle l\rangle}^\dag
    \wedge
    A_{\langle k\rangle}^\dag
  }{
    omitted; see PDF
  }
  
  and 

  \deqn{
    \left(
    A_{\langle k\rangle} B_{\langle l\rangle}
    \right)^\dag
    = 
    B_{\langle l\rangle}^\dag
    A_{\langle k\rangle}^\dag
  }{
    omitted; see PDF
  }

  See examples for package idiom.
  
}
\author{Robin K. S. Hankin}
\seealso{\code{\link{grade}},\code{\link{rev}}}
\references{C. Perwass 2009. \dQuote{Geometric algebra with applications
    in engineering}.  Springer.
}
\examples{

signature(2)

A <- rblade(g=3)
B <- rblade(g=4)

Conj(A \%^\% B) - Conj(B) \%^\% Conj(A)      # should be small
Conj(A * B) - Conj(B) * Conj(A)          # should be small


x1 <- rblade(d=9,g=2)
x2 <- rblade(d=9,g=2)
x3 <- rblade(d=9,g=2)
x4 <- rblade(d=9,g=2)

LHS <- Conj(x1 \%^\% x2 \%^\% x3 \%^\% x4)
RHS <- Conj(x4) \%^\% Conj(x3) \%^\% Conj(x2) \%^\% Conj(x1) 
Mod(LHS - RHS) # should be small


LHS <- Conj(x1 * x2 * x3 * x4)
RHS <- Conj(x4) * Conj(x3) * Conj(x2) * Conj(x1) 
Mod(LHS - RHS) # should be small


signature(0)

}