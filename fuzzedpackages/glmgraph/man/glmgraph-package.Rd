\name{glmgraph-package}
\alias{glmgraph-package}
\docType{package}
\title{Fit a GLM with a combination of sparse and smooth regularization}
\description{
    Fit a generalized linear model at grids of tuning parameter via penalized maximum likelihood.
    The regularization path is computed for a combination of sparse and smooth penalty 
    at two grids of values for the regularization parameter lambda1(Lasso or MCP penalty) and lambda2(Laplacian penalty).
    Fits linear, logistic regression models.
}
\details{
  \tabular{ll}{
    Package: \tab glmgraph\cr
    Type: \tab Package\cr
    Version: \tab 1.0-0\cr
    Date: \tab 2015-03-11\cr
    License: \tab GPL-2\cr}
  The algorithm accepts a design matrix \code{X}, a vector of responses \code{Y} and a Laplacian matrix \code{L}.
  Produces the regularization path over the grid of tuning parameter \code{lambda1} and \code{lambda2}. 
  It consists of the following main functions \cr
  \code{glmgraph}\cr
  \code{cv.glmgraph}\cr
  \code{plot.glmgraph}\cr
  \code{coef.glmgraph}\cr
  \code{predict.glmgraph}\cr
}
\references{Li Chen. Han Liu. Hongzhe Li. Jun Chen(2015) glmgraph: Graph-constrained Regularization for Sparse Generalized Linear Models.(Working paper)}
\author{Li Chen <li.chen@emory.edu>, Jun Chen <jun.chen2@mayo.edu>}
\keyword{models}
\keyword{regression}
\keyword{package}
\examples{
 set.seed(1234)
 library(glmgraph)
 n <- 100
 p1 <- 10
 p2 <- 90
 p <- p1+p2
 X <- matrix(rnorm(n*p), n,p)
 magnitude <- 1
 ## Construct Adjacency and Laplacian matrices
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
 Y <- eta+rnorm(n)
 obj <- glmgraph(X,Y,L,family="gaussian")
 plot(obj)
 betas <- coef(obj)
 betas <- coef(obj,lambda1=c(0.1,0.2))
 yhat <- predict(obj,X,type="response")
 cv.obj <- cv.glmgraph(X,Y,L)
 plot(cv.obj)
 beta.min <- coef(cv.obj)
 yhat.min <- predict(cv.obj,X)
 
}






