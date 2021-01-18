\name{APML0-package}
\alias{APML0-package}
\docType{package}
\title{
Augmented and Penalized Minimization Method L0
}
\description{
Fit linear, logistic and Cox models regularized with L0, lasso (L1), elastic-net (L1 and L2), or net (L1 and Laplacian) penalty, and their adaptive forms, such as adaptive lasso / elastic-net and net adjusting for signs of linked coefficients. It solves L0 penalty problem by simultaneously selecting regularization parameters and performing hard-thresholding (or selecting number of non-zeros). This augmented and penalized minimization method provides an approximation solution to the L0 penalty problem, but runs as fast as L1 regularization problem.

The package uses one-step coordinate descent algorithm and runs extremely fast by taking into account the sparsity structure of coefficients. It could deal with very high dimensional data.
}
\details{
\tabular{ll}{
Package: \tab APML0\cr
Type: \tab Package\cr
Version: \tab 0.10\cr
Date: \tab 2020-1-19\cr
License: \tab GPL (>= 2)\cr
}
Functions:
\code{\link{APML0}}, \code{\link{print.APML0}}\cr
}
\author{
Xiang Li, Shanghong Xie, Donglin Zeng and Yuanjia Wang\cr
Maintainer: Xiang Li <spiritcoke@gmail.com>
}
\references{Li, X., Xie, S., Zeng, D., Wang, Y. (2018).
\emph{Efficient l0-norm feature selection based on augmented and penalized minimization. Statistics in medicine, 37(3), 473-486.}\cr
\url{https://onlinelibrary.wiley.com/doi/full/10.1002/sim.7526}\cr
Boyd, S., Parikh, N., Chu, E., Peleato, B., Eckstein, J. (2011).
  \emph{Distributed optimization and statistical learning via the alternating direction method of multipliers. Foundations and Trends in Machine Learning, 3(1), 1-122.}\cr
  \url{http://dl.acm.org/citation.cfm?id=2185816}\cr
  Friedman, J., Hastie, T., Tibshirani, R. (2010).
  \emph{Regularization paths for generalized linear models via coordinate descent, Journal of Statistical Software, Vol. 33(1), 1.}\cr
  \url{http://www.jstatsoft.org/v33/i01/}\cr
}

\keyword{False positive control}
\keyword{Number of non-zeros}
\keyword{Regularization}
\keyword{Package}

\examples{
###  Linear model  ###
set.seed(1213)
N=100;p=30;p1=5
x=matrix(rnorm(N*p),N,p)
beta=rnorm(p1)
xb=x[,1:p1]\%*\%beta
y=rnorm(N,xb)

fiti=APML0(x,y,penalty="Lasso",nlambda=10) # Lasso
fiti2=APML0(x,y,penalty="Lasso",nlambda=10,nfolds=10) # Lasso
# attributes(fiti)



###  Logistic model  ###
set.seed(1213)
N=100;p=30;p1=5
x=matrix(rnorm(N*p),N,p)
beta=rnorm(p1)
xb=x[,1:p1]\%*\%beta
y=rbinom(n=N, size=1, prob=1.0/(1.0+exp(-xb)))

fiti=APML0(x,y,family="binomial",penalty="Lasso",nlambda=10) # Lasso
fiti2=APML0(x,y,family="binomial",penalty="Lasso",nlambda=10,nfolds=10) # Lasso
# attributes(fiti)



###  Cox model  ###
set.seed(1213)
N=100;p=30;p1=5
x=matrix(rnorm(N*p),N,p)
beta=rnorm(p1)
xb=x[,1:p1]\%*\%beta
ty=rexp(N, exp(xb))
td=rexp(N, 0.05)
tcens=ifelse(td<ty,1,0)  # censoring indicator
y=cbind(time=ty,status=1-tcens)

fiti=APML0(x,y,family="cox",penalty="Lasso",nlambda=10) # Lasso
fiti2=APML0(x,y,family="cox",penalty="Lasso",nlambda=10,nfolds=10) # Lasso
# attributes(fiti)
}


