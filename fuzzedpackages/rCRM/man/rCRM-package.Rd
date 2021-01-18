\name{rCRM-package}
\alias{rCRM-package}
\docType{package}
\title{
Regularized CRM
}
\description{
Fit a 2-parameter CRM model (O'Quigley and Shen 1996) regularized with L2 norm (Friedman et al. 2010) adjusted by the distance with the target DLT rate.

The package uses one-step coordinate descent algorithm and runs extremely fast.
}
\details{
\tabular{ll}{
Package: \tab rCRM\cr
Type: \tab Package\cr
Version: \tab 0.1\cr
Date: \tab 2018-11-06\cr
License: \tab GPL (>= 2)\cr
}
Functions:
\code{\link{rCRM}}\cr
}
\references{O'Quigley, J., Shen, L.Z. (1996).
\emph{Continual reassessment method: a likelihood approach. Biometrics, pp.673-684.}\cr
Friedman, J., Hastie, T. and Tibshirani, R. (2010).
\emph{Regularization paths for generalized linear models via coordinate descent, Journal of Statistical Software, Vol. 33(1), 1.}\cr
}
\author{
Xiang Li, Hong Tian, Kevin Liu, Pilar Lim\cr
Maintainer: Xiang Li <xli256@its.jnj.com>
}

\keyword{CRM}
\keyword{L2-norm}
\keyword{Regularization}
\keyword{Package}

\examples{

set.seed(1213)

dose0=c(1:6)
prob0=c(0.007, 0.086, 0.294, 0.545, 0.731, 0.841)

m=3; Y=NULL; X=NULL
for (i in 1:length(dose0)) {
  Y=c(Y, rbinom(m, size=1, prob=prob0[i]))
  X=c(X, rep(i, m))
}

fiti=rCRM(X, Y, dose0, tp=0.3, mlambda=10)
# attributes(fiti)


}


