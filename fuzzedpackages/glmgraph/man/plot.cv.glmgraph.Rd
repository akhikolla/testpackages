\name{plot.cv.glmgraph}
\alias{plot.cv.glmgraph}
\title{Plot the cross-validation curve produced by cv.glmgraph}
\description{Plots the cross-validation curve for the "cv.glmgraph"
  object, along with standard error bars.}
\usage{
\method{plot}{cv.glmgraph}(x,...) 
}
\arguments{
  \item{x}{A \code{"cv.glmgraph"} object.}
  \item{\dots}{Other graphical parameters to \code{plot}}
}
\references{Li Chen. Han Liu. Hongzhe Li. Jun Chen. (2015) glmgraph: Graph-constrained Regularization for Sparse Generalized Linear Models.(Working paper)}
\author{Li Chen <li.chen@emory.edu> , Jun Chen <chen.jun2@mayo.edu>}
\seealso{\code{glmgraph}, \code{cv.glmgraph}}
\examples{
 set.seed(1234)
 library(glmgraph)
 n <- 100
 p1 <- 10
 p2 <- 90
 p <- p1+p2
 X <- matrix(rnorm(n*p), n,p)
 magnitude <- 1
 A <- matrix(rep(0,p*p),p,p)
 A[1:p1,1:p1] <- 1
 A[(p1+1):p,(p1+1):p] <- 1
 diag(A) <- 0
 btrue <- c(rep(magnitude,p1),rep(0,p2))
 intercept <- 0
 eta <- intercept+X\%*\%btrue
 ### construct laplacian matrix from adjacency matrix
 diagL <- apply(A,1,sum)
 L <- -A
 diag(L) <- diagL
 ### gaussian
 Y <- eta+rnorm(n)
 cv.obj <- cv.glmgraph(X,Y,L)
 plot(cv.obj)
}
\keyword{models}
\keyword{regression}
