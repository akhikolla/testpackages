\name{print.cv.glmgraph}
\alias{print.cv.glmgraph}
\title{print a glmgraph object}
\description{
Print a summary of the cv.glmgraph solution path information during cross validation
}
\usage{
\method{print}{cv.glmgraph}(x, ...)
}
\arguments{
  \item{x}{fitted \code{cv.glmgraph} object}
  \item{\dots}{Other parameters to \code{print}}
}
\details{
 	The call prints the \code{cvmat} object from a fitted \code{cv.glmgraph} object.
 	The call also prints the chosen regularization parameters lambda1 and lambda2 along with best \code{cv.type}(minimum "mse" or "mae" if \code{family} is "gaussian"; maximum "auc" or minimum "deviance" if \code{family} is "binomial") after cross validation.
 }
\references{Li Chen. Han Liu. Hongzhe Li. Jun Chen. (2015) glmgraph: Graph-constrained Regularization for Sparse Generalized Linear Models.(Working paper)}
\author{Li Chen <li.chen@emory.edu> , Jun Chen <chen.jun2@mayo.edu>}
\seealso{\code{cv.glmgraph}}
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
 print(cv.obj)
}
\keyword{models}
\keyword{regression}
