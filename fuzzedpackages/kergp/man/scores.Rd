\name{scores}
\alias{scores}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Generic Function: Scores for a Covariance Kernel Object
}
\description{
  
  Generic function returning the scores for a covariance kernel object.

}
\usage{
scores(object, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{

    A covariance object.

  }
  \item{\dots}{
    
    Other arguments passed to methods.

  }
}
\details{

  Compute the derivatives \eqn{\partial_{\theta_k}\ell }{dl/d\theta[k]}
  for the (possibly concentrated) log-likelihood \eqn{\ell :=
  \log L}{l := log L} of a covariance object with parameter vector
  \eqn{\boldsymbol{\theta}}{\theta}.  The score for
  \eqn{\theta_k}{\theta[k]} is obtained as a matrix scalar product
  \deqn{
    \partial_{\theta_k} \ell
    = \textrm{trace}(\mathbf{W} \mathbf{D})
  }{dl/d\theta[k] = trace(W \%*\% D)}
  where \eqn{\mathbf{D} := \partial_{\theta_k} \mathbf{C}}{D = dC/d\theta[k]}
  and where  \eqn{\mathbf{W}}{W} is the matrix
  \eqn{
    \mathbf{W} := \mathbf{e}\mathbf{e}^\top - \mathbf{C}^{-1}
  }{e \%*\% t(e) - inv(C)}.
  The vector \eqn{\mathbf{e}}{e} is the vector of residuals
  and the matrix \eqn{\mathbf{C}}{C}
  is the covariance computed for the design \eqn{\mathbf{X}}{X}.
  
}
\value{

  A numeric vector of length \code{npar(object)} containing the scores.

}
%%\references{}
%%\author{}
\note{

  The scores can be efficiently computed when the matrix
  \eqn{\mathbf{W}}{W} has already been pre-computed.

}
%%\seealso{}
%%\examples{}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%%\keyword{ ~kwd1 }
%%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
