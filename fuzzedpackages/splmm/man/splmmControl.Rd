\name{splmmControl}
\alias{splmmControl}
\title{Options for the 'splmm' Algorithm}
\description{Definition of various kinds of options in the
  algorithm.}
\usage{splmmControl(tol=10^(-4),trace=1,maxIter=1000,maxArmijo=20,number=5,a_init=1,
           delta=0.1,rho=0.001,gamma=0,lower=10^(-6),upper=10^8,seed=532,VarInt=c(0,10),
           CovInt=c(-5,5),thres=10^(-4))}
\arguments{
  \item{tol}{convergence tolerance}
  \item{trace}{integer. 1 prints no output, 2 prints warnings, 3 prints
    the current function values and warnings (not recommended)}
  \item{maxIter}{maximum number of (outer) iterations}
  \item{maxArmijo}{maximum number of steps to be chosen in the Armijo Rule. If the maximum is reached, the algorithm continues with optimizing the next coordinate.}
  \item{number}{integer. Determines the active set algorithm. The zero
    fixed-effects coefficients are only updated each number
    iteration. It may be that a smaller number increases the speed of
    the algorithm. Use \eqn{0 \le number \le 5}.}
  \item{a_init}{\eqn{\alpha_{init}} in the Armijo step. See Schelldorfer et. al. (2010).}
  \item{delta}{\eqn{\delta} in the Armijo step. See Schelldorfer et. al. (2010)}
  \item{rho}{\eqn{\rho} in the Armijo step. See Schelldorfer et. al. (2010)}
  \item{gamma}{\eqn{\gamma} in the Armijo step. See Schelldorfer et. al. (2010)}
  \item{lower}{lower bound for the Hessian}
  \item{upper}{upper bound for the Hessian}
  \item{seed}{set.seed for calculating the starting value, which
    performs a 10-fold cross-validation.}
  \item{VarInt}{Only for opt="optimize". The interval for the variance
    parameters used in "optimize". See help("optimize")}
  \item{CovInt}{Only for opt="optimize". The interval for the covariance
    parameters used in "optimize". See help("optimize")}
  \item{thres}{If a variance or covariance parameter has smaller absolute value than thres, the parameter is set to exactly zero.}
}
\details{For the Armijo step parameters, see Bertsekas (2003)}



% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{misc}
