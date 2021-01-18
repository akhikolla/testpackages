\name{coef.cv.glmgraph}
\alias{coef.cv.glmgraph}
\title{Retrieve coefficients from a fitted "cv.glmgraph" object.}
\description{
Retrieve coefficients from a fitted "cv.glmgraph" object based on the chosen regularization parameters from cross validation.}
\usage{
\method{coef}{cv.glmgraph}(object,s,...)
}
\arguments{
  \item{object}{Fitted \code{"cv.glmgraph"} model object.}
  \item{s}{Either \code{"lambda1.min"} or \code{"lambda1.1se"}. If \code{"lambda1.min"} is used, coefficients of best cross validation criteria(minimum "mse" or "mae" if \code{family} is "gaussian"; maximum "auc" or minimum "deviance" if \code{family} is "binomial") are returned. Otherwise, coefficients based on one-standard error rule are returned. The default value is \code{"lambda1.min"}.}
  \item{\dots}{Other parameters to \code{coef}}
 }
\references{Li Chen. Han Liu. Hongzhe Li. Jun Chen. (2015) glmgraph: Graph-constrained Regularization for Sparse Generalized Linear Models.(Working paper)}
\author{Li Chen <li.chen@emory.edu> , Jun Chen <chen.jun2@mayo.edu>}
\seealso{\code{predict.cv.glmgraph},\code{cv.glmgraph}}
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
 diagL <- apply(A,1,sum)
 L <- -A
 diag(L) <- diagL
 btrue <- c(rep(magnitude,p1),rep(0,p2))
 intercept <- 0
 eta <- intercept+X\%*\%btrue 
 ### gaussian
 Y <- eta+rnorm(n)
 cv.obj <- cv.glmgraph(X,Y,L)
 beta.min <- coef(cv.obj)
}
\keyword{models}
\keyword{regression}
