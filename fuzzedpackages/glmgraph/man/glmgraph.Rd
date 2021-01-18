\name{glmgraph}
\alias{glmgraph}
\title{Fit a GLM with a combination of sparse and smooth regularization}
\description{
    Fit a generalized linear model at grids of tuning parameter via penalized maximum likelihood.
    The regularization path is computed for a combination of sparse and smooth penalty 
    at two grids of values for the regularization parameter lambda1(Lasso or MCP penalty) and lambda2(Laplacian penalty).
    Fits linear, logistic regression models.
}

\usage{
glmgraph(X, Y, L, family=c("gaussian","binomial"), penalty=c("MCP","lasso") ,
mcpapproach=c("mmcd", "adaptive", "original"),gamma=8,
lambda1,nlambda1=100,lambda2=c(0, 0.01 * 2^(0:7)),eps=1e-3,max.iter=2000,
dfmax=round(ncol(X)/2),penalty.factor=rep(1,ncol(X)),standardize=TRUE,warn=FALSE,...)
}

\arguments{
  \item{X}{Input matrix; each row is an observation vector.}
  \item{Y}{Response vector. Quantitative for \code{family="gaussian"} or binary(0/1) for \code{family="binomial"}.}
  \item{family}{Either "gaussian", "binomial", depending
    on the response.}
  \item{L}{User-specified Laplacian matrix.}
  \item{penalty}{The sparse penalty to be applied to the model.  Either "MCP"
    (the default), or "lasso".}
   \item{mcpapproach}{For \code{family="binomial"}, three optional algorithms are provided when \code{penalty} is set to MCP: 
    "mmcd"(Majorization minimization by coordinate descent);
    "adaptive"(Adaptive rescaling) and "original"(without any adjustment).
    For \code{family="gaussian"}, the option could only be "original".
    }
  \item{gamma}{The tuning parameter of the MCP penalty.The default value is 8.}
  \item{nlambda1}{The number of \code{lambda1} values.  Default is 100.}
  \item{lambda1}{A user-specified sequence of \code{lambda1} values. Typical usage
    is to have the program compute its own \code{lambda1} sequence based on
    \code{nlambda1} and \code{lambda1.min.ratio}. Supplying a value of \code{lambda1} overrides this.
  	By default, a sequence of values of length \code{nlambda1} is computed, equally
    spaced on the log scale.}
  \item{lambda2}{A user-specified sequence of \code{lambda2} values. The default value
  	are 0 and 0.01*2^(0:7). The selection of \code{lambda2} depends on the data and should be adapted in some cases. 
  	A good suggestion is to try a few \code{lambda2} and plot the results. } 
  \item{eps}{ Convergence threshold for coordinate descent. Each inner
  	coordinate-descent loop continues until the relative change in the
  	objective function is less than \code{eps1}. 
    Default is \code{1e-3}.}
  \item{max.iter}{Maximum number of passes over the data for all \code{lambda1}
  values. Default is 2000.}
  \item{dfmax}{Limit the maximum number of variables in the
    model. Useful for very large \code{p}.
    Default value equals to half of \code{p}.}
  \item{penalty.factor}{A multiplicative factor for the penalty applied
    to each coefficient. If supplied, \code{penalty.factor} must be a
    numeric vector of length equal to the number of columns of
    \code{X}.  The purpose of \code{penalty.factor} is to apply
    differential penalization if some coefficients are thought to be
    more likely than others to be in the model.  In particular,
    \code{penalty.factor} can be 0, in which case the coefficient is
    always in the model without shrinkage.}
  \item{standardize}{
  	Logical flag for x variable standardization, prior to
    fitting the model sequence. The coefficients are always returned on
    the original scale. Default is \code{standardize=TRUE}.
 	If variables are in the same units already, you might not wish to standardize.
    }
  \item{warn}{Return warning messages for failures to converge and model selection issues. Default is FALSE.}
  \item{\dots}{Other parameters to \code{glmgraph}}
}

\value{
  An object \code{"glmgraph"} containing:
  \item{betas}{A list of fitted coefficients. The number of rows for each matrix is equal to the number of coefficients, and the number of columns is smaller or equal to \code{nlambda1}.}
   \item{lambda1s}{A list of vector. Each vector is a sequence of used \code{lambda1} for each used \code{lambda2}.}
   \item{lambda2}{A sequence of \code{lambda2} actually used.}
   \item{loglik}{A list of log likelihood for each value of \code{lambda1} and \code{lambda2}.}
   \item{df}{A list of the number of nonzero values for each value of \code{lambda1} and \code{lambda2}.}
}

\references{Li Chen. Han Liu. Hongzhe Li. Jun Chen. (2015) Graph-constrained Regularization for Sparse Generalized Linear Models.(Working paper)}
\author{Li Chen <li.chen@emory.edu>, Jun Chen <chen.jun2@mayo.edu>}
\seealso{\code{plot.glmgraph}, \code{cv.glmgraph}}
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
 obj <- glmgraph(X,Y,L,family="gaussian")
 plot(obj)
 ### binomial
 Y <- rbinom(n,1,prob=1/(1+exp(-eta)))
 obj <- glmgraph(X,Y,L,family="binomial")
 plot(obj) 
 
}
\keyword{models}
\keyword{regression}
