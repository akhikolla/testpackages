\name{plot.glmgraph}
\alias{plot.glmgraph}
\title{Plot coefficients from a "glmgraph" object}
\description{Plot solution path for a fitted
  \code{"glmgraph"} object.}
\usage{
\method{plot}{glmgraph}(x,...)
}
\arguments{
  \item{x}{Fitted \code{"glmgraph"} model.}
  \item{\dots}{Other graphical parameters to \code{plot}}
}
\references{Li Chen. Han Liu. Hongzhe Li. Jun Chen. (2015) glmgraph: Graph-constrained Regularization for Sparse Generalized Linear Models.(Working paper)}
\author{Li Chen <li.chen@emory.edu> , Jun Chen <chen.jun2@mayo.edu>}
\seealso{\code{glmgraph}}
\examples{
 set.seed(1234)
 library(glmgraph)
 n <- 100
 p1 <- 10
 p2 <- 90
 p <- p1+p2
 X <- matrix(rnorm(n*p), n,p)
 magnitude <- 1
 ### construct laplacian matrix from adjacency matrix
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
 plot(obj)
 ### binomial
 Y <- rbinom(n,1,prob=1/(1+exp(-eta)))
 obj <- glmgraph(X,Y,L,family="binomial")
 plot(obj) 

}
\keyword{models}
\keyword{regression}
