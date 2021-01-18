\name{deriv}
\alias{deriv}
\alias{deriv.spray}
\alias{aderiv}
\title{
Partial differentiation of spray objects
}
\description{\loadmathjax
  Partial differentiation of spray objects interpreted as multivariate
  polynomials
}
\usage{
\method{deriv}{spray}(expr, i , derivative = 1, ...)
aderiv(S,orders)
}
\arguments{
  \item{expr}{A spray object, interpreted as a multivariate polynomial}
  \item{i}{Dimension to differentiate with respect to}
  \item{derivative}{How many times to differentiate}
  \item{...}{Further arguments, currently ignored}
  \item{S}{spray object}
  \item{orders}{The orders of the differentials}
}
\details{\loadmathjax

Function \code{deriv.spray()} is the method for generic \code{spray()};
if \code{S} is a spray object, then \code{spray(S,i,n)} returns
\mjeqn{\partial^n S/\partial x_i^n =
S^{\left(x_i,\ldots,x_i\right)}}{d^n S/dx_i^n = S^(x_i...x_i)}.

Function \code{aderiv()} is the generalized derivative; if \code{S} is a
spray of arity 3, then \code{aderiv(S,c(i,j,k))} returns
\mjeqn{\frac{\partial^{i+j+k} S}{\partial x_1^i\partial x_2^j\partial
x_3^k}}{d^(i+j+k)S/dx_1^i dx_2^j dx_3^k}.

}
\value{
Both functions return a spray object.
}
\author{
Robin K. S. Hankin
}
\seealso{\code{\link{asum}}}

\examples{


S <- spray(matrix(sample(-2:2,15,replace=TRUE),ncol=3),addrepeats=TRUE)

deriv(S,1)
deriv(S,2,2)

# differentiation is invariant under order:
aderiv(S,1:3) == deriv(deriv(deriv(S,1,1),2,2),3,3)

# Leibniz's rule:
S1 <- spray(matrix(sample(0:3,replace=TRUE,21),ncol=3),sample(7),addrepeats=TRUE)
S2 <- spray(matrix(sample(0:3,replace=TRUE,15),ncol=3),sample(5),addrepeats=TRUE)

S1*deriv(S2,1) + deriv(S1,1)*S2 == deriv(S1*S2,1)

# Generalized Leibniz:
aderiv(S1*S2,c(1,1,0)) == (
aderiv(S1,c(0,0,0))*aderiv(S2,c(1,1,0)) +
aderiv(S1,c(0,1,0))*aderiv(S2,c(1,0,0)) +
aderiv(S1,c(1,0,0))*aderiv(S2,c(0,1,0)) +
aderiv(S1,c(1,1,0))*aderiv(S2,c(0,0,0)) 
)



}

\keyword{mathsymbol}

