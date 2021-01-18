\name{coef.glmgraph}
\alias{coef.glmgraph}
\title{Retrieve coefficients from a fitted "glmgraph" object.}
\description{
Retrieve coefficients from a fitted "glmgraph" object, depending on the user-specified regularization parameters.
}
\usage{
\method{coef}{glmgraph}(object,lambda1,lambda2,...)
}
\arguments{
  \item{object}{Fitted \code{"glmgraph"} model object.}
  \item{lambda1}{Values of the regularization parameter \code{lambda1} at
    which retrieval of coefficients are requested. For values of \code{lambda1} not in
    the sequence of fitted models, linear interpolation is used. 
    However, \code{lambda1} should be within the range of \code{lambda1} used to fit \code{glmgraph} object.}
  \item{lambda2}{The user-specified regularization \code{lambda2} should be exactly 
  	subset of the \code{lambda2} used to fit \code{glmgraph} object. Linear interpolation is not used.}
  \item{\dots}{Other parameters to \code{coef}}
}
\value{The object returned depends on type.}
\details{
	If \code{lambda1} and \code{lambda2} are missing, all coefficients of fitted \code{glmgraph} object will be returned.
    If only \code{lambda1} is missing, then coefficients of specified \code{lambda2} will be returned.
}
\references{Li Chen. Han Liu. Hongzhe Li. Jun Chen. (2015) glmgraph: Graph-constrained Regularization for Sparse Generalized Linear Models.(Working paper)}
\author{Li Chen <li.chen@emory.edu>, Jun Chen <chen.jun2@mayo.edu>}
\seealso{\code{predict.glmgraph},\code{glmgraph}}
\examples{
 set.seed(1234)
 library(glmgraph)
 n <- 100
 p1 <- 10
 p2 <- 90
 p <- p1+p2
 X <- matrix(rnorm(n*p), n,p)
 magnitude <- 1
 ## construct laplacian matrix from adjacency matrix
 A <- matrix(rep(0,p*p),p,p)
 A[1:p1,1:p1] <- 1
 A[(p1+1):p,(p1+1):p] <- 1
 diag(A) <- 0
 btrue <- c(rep(magnitude,p1),rep(0,p2))
 intercept <- 0
 eta <- intercept+X\%*\%btrue
 diagL <- apply(A,1,sum)
 L <- -A
 diag(L) <- diagL
 ### gaussian
 Y <- eta+rnorm(n)
 obj <- glmgraph(X,Y,L)
 coefs <- coef(obj)
 coefs <- coef(obj,lambda2=0.01)
 coefs <- coef(obj,lambda1=c(0.11,0.12))
 coefs <- coef(obj,lambda1=c(0.11,0.12),lambda2=0.01)
}
\keyword{models}
\keyword{regression}
