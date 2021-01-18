\name{ADMMnet-package}
\alias{ADMMnet-package}
\docType{package}
\title{
Regularized Model with Selecting the Number of Non-Zeros
}
\description{
This package fits linear and cox models regularized with net (L1 and Laplacian), elastic-net (L1 and L2) or lasso (L1) penalty, and their adaptive forms, such as adaptive lasso and net adjusting for signs of linked coefficients. In addition, it treats the number of non-zero coefficients as another tuning parameter and simultaneously selects with the regularization parameter \code{lambda}. This is motivated by formulating L0 variable selection in ADMM form. By selecting the regularization parameter and the number of non-zeros, it shows significant improvement over the commonly used regularized methods, which depend on the regularization parameter only.

The package uses one-step coordinate descent algorithm and runs extremely fast by taking into account the sparsity structure of coefficients.
}
\details{
\tabular{ll}{
Package: \tab ADMMnet\cr
Type: \tab Package\cr
Version: \tab 0.1\cr
Date: \tab 2015-12-10\cr
License: \tab GPL (>= 2)\cr
}
Functions:
\code{\link{ADMMnet}}, \code{\link{print.ADMMnet}}\cr
}
\author{
Xiang Li, Shanghong Xie, Donglin Zeng and Yuanjia Wang\cr
Maintainer: Xiang Li <xl2473@cumc.columbia.edu>
}
\references{Boyd, S., Parikh, N., Chu, E., Peleato, B., & Eckstein, J. (2011).
  \emph{Distributed optimization and statistical learning via the alternating direction method of multipliers. Foundations and Trends in Machine Learning, 3(1), 1-122.}\cr
  \url{http://dl.acm.org/citation.cfm?id=2185816}\cr
  Friedman, J., Hastie, T. and Tibshirani, R. (2010).
  \emph{Regularization paths for generalized linear models via coordinate descent, Journal of Statistical Software, Vol. 33(1), 1.}\cr
  \url{http://www.jstatsoft.org/v33/i01/}\cr
  Li, C., and Li, H. (2010).
  \emph{Variable selection and regression analysis for graph-structured covariates with an application to genomics. The annals of applied statistics, 4(3), 1498.}\cr
  \url{http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3423227/}\cr
   Sun, H., Lin, W., Feng, R., and Li, H. (2014)
  \emph{Network-regularized high-dimensional cox regression for analysis of genomic data, Statistica Sinica.}\cr
  \url{http://www3.stat.sinica.edu.tw/statistica/j24n3/j24n319/j24n319.html}\cr
}

\keyword{Number of non-zeros}
\keyword{Regularization}
\keyword{Package}

\examples{
###  Linear model  ###
set.seed(1213)
N=100;p=30;p1=5
x=matrix(rnorm(N*p),N,p)
beta=rnorm(p1)
xb=x[,1:p1]%*%beta
y=rnorm(N,xb)

fiti=ADMMnet(x,y,penalty="Lasso",nlambda=10,nfolds=10) # Lasso
# attributes(fiti)


###  Cox model  ###
set.seed(1213)
N=100;p=30;p1=5
x=matrix(rnorm(N*p),N,p)
beta=rnorm(p1)
xb=x[,1:p1]%*%beta
ty=rexp(N,exp(xb))
tcens=rbinom(n=N,prob=.3,size=1)  # censoring indicator
y=cbind(time=ty,status=1-tcens)

fiti=ADMMnet(x,y,family="cox",penalty="Lasso",nlambda=10,nfolds=10) # Lasso
# attributes(fiti)
}


