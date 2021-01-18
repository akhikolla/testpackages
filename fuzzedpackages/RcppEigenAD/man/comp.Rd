\docType{methods}
\name{\%.\%}
\alias{\%.\%}
\title{Compose two functions created using either \link{sourceCppAD} or \code{\%.\%} itself }
\usage{
f \%.\% g
}
\arguments{
\item{f}{Function to be composed with \code{g}}

\item{g}{Function to be composed with \code{f}}
}
\value{
A function which computes the composition \code{f} and \code{g}
}
\description{
Returns a function with a matrix \code{X} as input which computes the value of \eqn{(f \circ g)(x) = f(g(x))}. 
Note that the order of the composition is such that \code{g} is applied to \code{X} first. \code{f} is then applied
to the result of \code{g}. The returned function is compatible with both \link{J} and \link{H}, and if \link{J} and \link{H} are applied to
a function produced by composition, the resulting Jacobian or Hessian matrices are constructed
from the Jacobians and Hessians of \code{f} and \code{g} using a
combinatorical form of Faa di Bruno's formula (Hardy 2006). The functions \code{f} and \code{g} must both be functions of a single argument.
Note that a function of multiple arguments can be Curried into a function with a single argument. The R package functional provides
the method \code{\link[functional]{Curry}} which is convenient for this purpose.
}
\examples{
\donttest{
library(RcppEigenAD)
# define a function to calculate the eigen vectors of a matrix
f<-sourceCppAD('
ADmat f(const ADmat& X)
{
  Eigen::EigenSolver<ADmat> es(X);
  return es.pseudoEigenvectors();

}
')
# define function to calculate the inverse of a matrix
g<-sourceCppAD('
ADmat g(const ADmat& X)
{
  return X.inverse(); 
}
')
# compose f and g to produce a functions to calculate the eigenvectors of the inverse of a matrix
h<-f\%.\%g  #  h = f o g
x<-matrix(c(1,2,3,4),2,2)
x<-x\%*\%t(x) # positive definite matrix
J(h)(x) # Jacobian of h = f o g
H(h)(x) # Stacked Hessians of h = f o g
# redefine h as a function to directly calculate the eigenvectors of the inverse of a matrix
h<-sourceCppAD('
ADmat h(const ADmat& X)
{
  Eigen::EigenSolver<ADmat> es(X.inverse());
  return es.pseudoEigenvectors();
}
')
# calculate the Jacobian and Hessian of h to compare with previous result
J(h)(x) # Jacobian of h = f o g
H(h)(x) # Stacked Hessians of h = f o g
}
}
\references{
  \insertRef{Hardy06}{RcppEigenAD}
  
  \insertRef{ma-tsoy-wo-faa-di-bruno}{RcppEigenAD}
}
