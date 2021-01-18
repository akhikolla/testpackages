\name{cv.glmgraph}
\alias{cv.glmgraph}
\title{Cross-validation for glmgraph}
\description{Performs k-fold cross validation for glmgraph}
\usage{
cv.glmgraph(X,Y,L,...,type.measure=c("mse","mae","deviance","auc"),nfolds=5,trace=TRUE)
}
\arguments{
  \item{X}{\code{X} matrix as in \code{glmgraph}.}
  \item{Y}{Response \code{Y} as in \code{glmgraph}.}
  \item{L}{User-specified Laplacian matrix \code{L} as in \code{glmgraph}.}
  \item{...}{Additional arguments as in \code{glmgraph}.}
  \item{type.measure}{if \code{family} is "gaussian", the \code{type.measure} option is "mse"(mean squared error) or "mae"(mean absolute error);
  if \code{family} is "binomial", the \code{type.measure} option is "deviance" or "auc"(area under the curve). The default is "mse".}
  \item{nfolds}{The number of cross-validation folds. Default is 5.}
  \item{trace}{Print out the cross validation steps if \code{trace} is specified TRUE.}
 }
 \details{The function runs \code{glmgraph} \code{nfolds}+1 times; the
  first to get the \code{lambda1} and \code{lambda2} sequence, and then the remainder to
  compute the fit with each of the folds omitted. The error is
  accumulated, and the average error and standard deviation over the
  folds is computed. 
  Note also that the results of \code{cv.glmgraph} are
  random, since the folds are selected at random. Users can reduce this
  randomness by running \code{cv.glmgraph} many times, and averaging the
  error curves. 
}

\value{	
  An object \code{"cv.glmgraph"} containing:
  \item{obj}{The fitted \code{glmgraph} object for the whole data.}
  \item{cvmat}{
  	A data frame summarized cross validation results, which could be obtained by \code{print} function.
  	It has \code{lambda2},\code{lambda1.min},\code{cvmin},\code{semin},\code{lambda1.1se} as columns. 
  	Each row represents that for this \code{lambda2}, \code{lambda1} with best \code{type.measure} \code{cvmin} is chosen and reported as \code{lambda1.min}. 
  	If one standard error rule is applied, \code{lambda1.1se} and its corresponding best \code{type.measure} value \code{semin} is reported.
  }
  \item{cvm}{The mean cross-validated \code{type.measure} value. A list of vector contains \code{type.measure}.
  	Each element of the list is a vector that is \code{type.measure} value for one \code{lambda2} 
  	across all \code{lambda1} sequence averaged across K-fold.}
  \item{cvsd}{The estimate of standard error of \code{cvm}.}
  \item{cvmin}{Best cross-validation \code{type.measure} value across all combination of \code{lambda1} and \code{lambda2}. 
  It is minimum "mse" or "mae" if \code{family} is "gaussian"; it is the maximum "auc" or minimum "deviance" if \code{family} is "binomial".}  
  \item{cv.1se}{Simliar to \code{cvmin} except one standard error rule is applied.}
  \item{lambda1.min}{Coupled with \code{lambda2.min} is the optimal regularization parameter selection.}
  \item{lambda2.min}{Coupled with \code{lambda1.min} is the optimal regularization parameter selection.}
  \item{lambda1.1se}{Coupled with \code{lambda2.1se} is the optimal regularization parameter selection if one standard error rule is applied.}
  \item{lambda2.1se}{Coupled with \code{lambda1.1se} is the optimal regularization parameter selection if one standard error rule is applied. }
  \item{beta.min}{Estimated beta with best \code{type.measure} value with the regularization parameter of \code{lambda1.min} and \code{lambda2.min}.}
  \item{beta.1se}{Estimated beta with best \code{type.measure} value with the regularization parameter of \code{lambda1.1se} and \code{lambda2.1se}.}
  }
\references{Li Chen. Han Liu. Hongzhe Li. Jun Chen(2015) glmgraph: Graph-constrained Regularization for Sparse Generalized Linear Models.(Working paper)}
\author{Li Chen <li.chen@emory.edu> , Jun Chen <chen.jun2@mayo.edu>}
\seealso{\code{glmgraph},\code{coef.cv.glmgraph},\code{predict.cv.glmgraph}}
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
 cv.obj <- cv.glmgraph(X,Y,L,penalty="lasso",lambda2=c(0,1.28))
 beta.min <- coef(cv.obj)
 print(cv.obj)
 ### binomial
 Y <- rbinom(n,1,prob=1/(1+exp(-eta)))
 cv.obj <- cv.glmgraph(X,Y,L,family="binomial",lambda2=c(0,1.28),penalty="lasso",type.measure="auc")
 beta.min <- coef(cv.obj)
 print(cv.obj)
}


