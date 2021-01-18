\name{predict.cv.glmgraph}
\alias{predict.cv.glmgraph}
\title{make prediction from a fitted "cv.glmgraph" object.}
\description{
This function makes predictions from a cross-validated glmgraph model,
using the stored \code{"cv.glmgraph"} object, and the optimal value
chosen for \code{lambda1} and \code{lambda2}.
}
\usage{
\method{predict}{cv.glmgraph}(object,X,s=c("lambda1.min","lambda1.1se"),
type=c("response", "coefficients","class", "nzeros","link"),...)
}
\arguments{
  \item{object}{Fitted \code{"cv.glmgraph"} model object.}
  \item{X}{Matrix at which predictions are to be made.}
  \item{s}{Either \code{"lambda1.min"} or \code{"lambda1.1se"}.If \code{"lambda1.min"} is used, prediction based on coefficient of best cross validation criteria(minimum "mse" or "mae" if \code{family} is "gaussian"; maximum "auc" or minimum "deviance" if \code{family} is "binomial") are returned. Otherwise, predictficients based on one-standard error rule are returned. The default value is \code{"lambda1.min"}.}
    \item{type}{Type of prediction: \code{"link"} returns the
    linear predictors; \code{"response"} gives the fitted values;
    \code{"class"} returns the binomial outcome with the highest
    probability; \code{"coefficients"} returns the coefficients;
    \code{"nzeros"} returns a list containing the indices and names of the
    nonzero variables at each combination of \code{lambda1} and \code{lambda2}.
	}
	\item{\dots}{Other parameters to \code{predict}}
 }

\references{Li Chen. Han Liu. Hongzhe Li. Jun Chen. (2015) Graph-constrained Regularization for Sparse Generalized Linear Models.(Working paper)}
\author{Li Chen <li.chen@emory.edu> , Jun Chen <chen.jun2@emory.edu>}
\seealso{\code{cv.glmgraph},\code{coef.cv.glmgraph}}
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
 cv.obj <- cv.glmgraph(X,Y,L)
 beta.min <- predict(cv.obj,X,type="coefficients")
}
\keyword{models}
\keyword{regression}
