\name{predict.glmgraph}
\alias{predict.glmgraph}
\alias{predict.glmgraph}
\title{Model predictions based on a fitted "glmgraph" object.}
\description{
Similar to other predict methods, this function returns predictions from
a fitted \code{"glmgraph"} object.
}
\usage{
\method{predict}{glmgraph}(object, X, type=c("response", "coefficients",
"class", "nzeros","link"), lambda1, lambda2,...)
}
\arguments{
  \item{object}{Fitted \code{"glmgraph"} model object.}
  \item{X}{Matrix of values at which predictions are to be made.}
  \item{lambda1}{Values of the regularization parameter \code{lambda1} at
    which predictions are requested.  For values of \code{lambda1} not in
    the sequence of fitted models, linear interpolation is used.}
    \item{lambda2}{ Values of the regularization parameter \code{lambda1} at
    which predictions are requested. Specified \code{lambda2} should be the subet
    of \code{lambda2} used to fit \code{glmgraph} object.}
    
    \item{type}{Type of prediction: \code{"link"} returns the
    linear predictors; \code{"response"} gives the fitted values;
    \code{"class"} returns the binomial outcome with the highest
    probability; \code{"coefficients"} returns the coefficients;
    \code{"nzeros"} returns a list containing the indices and names of the
    nonzero variables at each combination of \code{lambda1} and \code{lambda2}.
    }
    \item{\dots}{Other parameters to \code{predict}}
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
 res <- predict(obj, X, type="link", lambda1=0.05,lambda2=0.01)
 res <- predict(obj, X, type="response", lambda1=0.05,lambda2=0.01)
 res <- predict(obj,X,type="nzeros",lambda1=0.05,lambda2=0.01)
 ### binomial
 Y <- rbinom(n,1,prob=1/(1+exp(-eta)))
 obj <- glmgraph(X,Y,L,family="binomial")
 res <- predict(obj,X,type="class",lambda1=c(0.05,0.06),lambda2=c(0.02,0.16,0.32))
}
\keyword{models}
\keyword{regression}
