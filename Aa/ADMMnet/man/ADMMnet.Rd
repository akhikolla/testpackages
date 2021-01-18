\name{ADMMnet}
\alias{ADMMnet}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Fit a Model with Various Regularization Forms
}
\description{
Fit a linear or cox model regularized with net (L1 and Laplacian), elastic-net (L1 and L2) or lasso (L1) penalty, and their adaptive forms, such as adaptive lasso and net adjusting for signs of linked coefficients. In addition, it treats the number of non-zero coefficients as another tuning parameter and simultaneously selects with the regularization parameter \code{lambda}. The package uses one-step coordinate descent algorithm and runs extremely fast by taking into account the sparsity structure of coefficients.
}
\usage{
ADMMnet(x, y, family = c("gaussian", "cox"), penalty = c("Lasso", "Enet", "Net"),
  Omega = NULL, alpha = 1.0, lambda = NULL, nlambda = 50, rlambda = NULL,
  nfolds = 1, foldid = NULL, inzero = TRUE, adaptive = c(FALSE, TRUE), aini = NULL,
  isd = FALSE, keep.beta = FALSE, ifast = TRUE, thresh = 1e-07, maxit = 1e+05)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{input matrix. Each row is an observation vector.
}
  \item{y}{response variable. For \code{family = "gaussian"},
  \code{y} is a continuous vector. For \code{family = "cox"}, \code{y} is a two-column matrix with columns named `time' and `status'. `status' is a binary variable, with `1' indicating event, and `0' indicating right censored.
}
  \item{family}{type of outcome. Can be "gaussian" or "cox".
}
  \item{penalty}{penalty type. Can choose \code{"Net"}, \code{"Enet"} (elastic net) and \code{"Lasso"}. For \code{"Net"}, need to specify \code{Omega}; otherwises, \code{"Enet"} is performed. For \code{penalty = "Net"}, the penalty is defined as \deqn{\lambda*{\alpha*||\beta||_1+(1-\alpha)/2*(\beta^{T}L\beta)},}
where \eqn{L} is a Laplacian matrix calculated from \code{Omega}.
}
  \item{Omega}{correlation/adjacency matrix with zero diagonal, used for \code{penalty = "Net"} to calculate Laplacian matrix.
}
  \item{alpha}{ratio between L1 and Laplacian for \code{"Net"}, or between L1 and L2 for \code{"Enet"}. Default is \code{alpha = 1.0}, i.e. lasso.
}
  \item{lambda}{a user supplied decreasing sequence. If \code{lambda = NULL}, a sequency of \code{lambda} is generated based on \code{nlambda} and \code{rlambda}. Supplying a value of \code{lambda} overrides this.
}
  \item{nlambda}{number of \code{lambda} values. Default is 50.
}
  \item{rlambda}{fraction of \code{lambda.max} to determine the smallest value for \code{lambda}. The default is \code{rlambda = 0.0001} when the number of observations is larger than or equal to the number of variables; otherwise, \code{rlambda = 0.01}.
}
  \item{nfolds}{number of folds. With \code{nfolds = 1} and \code{foldid = NULL} by default, cross-validation is not performed. For cross-validation, smallest value allowable is \code{nfolds = 3}. Specifying \code{foldid} overrisdes \code{nfolds}.
}
  \item{foldid}{an optional vector of values between 1 and \code{nfolds} specifying which fold each observation is in.
}
  \item{inzero}{logical flag for simultaneously selecting the number of non-zero coefficients with \code{lambda}. Default is \code{inzero = TRUE}.
}
  \item{adaptive}{logical flags for adaptive version. Default is \code{adaptive = c(FALSE, TRUE)}. The first element is for adaptive on \eqn{\beta} in L1 and the second for adjusting for signs of linked coefficients in Laplacian matrix.
}
  \item{aini}{a user supplied initial estimate of \eqn{\beta}. It is a list including \code{wbeta} for adaptive L1 and \code{sgn} for adjusting Laplacian matrix. \code{wbeta} is the absolute value of inverse initial estimates. If \code{aini = NULL} but \code{adaptive} is required, \code{aini} is generated from regularized model with \code{penatly = "Enet"} and \code{alpha = 0.0}, i.e. a ridge regression.
}
  \item{isd}{logical flag for outputing standardized coefficients. \code{x} is always standardized prior to fitting the model. Default is \code{isd = FALSE}, returning \eqn{\beta} on the original scale.
}
  \item{keep.beta}{logical flag for returning estimates for all \code{lambda} values. For \code{keep.beta = FALSE}, only return the estimate with the minimum cross-validation value.
}
  \item{ifast}{logical flag for efficient calculation of risk set updates for \code{family = "cox"}. Default is \code{ifast = TRUE}.
}
  \item{thresh}{convergence threshold for coordinate descent. Default value is \code{1E-7}.
}
  \item{maxit}{Maximum number of iterations for coordinate descent. Default is \code{10^5}.
}
}
\details{
  One-step coordinate descent algorithm is applied for each \code{lambda}. For \code{family = "cox"}, \code{ifast = TRUE} adopts an efficient way to update risk set and sometimes the algorithm ends before all \code{nlambda} values of \code{lambda} have been evaluated. To evaluate small values of \code{lambda}, use \code{ifast = FALSE}. The two methods only affect the efficiency of algorithm, not the estimates.
  
  \code{x} is always standardized prior to fitting the model and the estimate is returned on the original scale. For \code{family = "gaussian"}, y is centered by removing its mean, so there is no intercept output.
  
  Cross-validation is used for tuning parameters. For \code{inzero =  TRUE}, we further select the number of non-zero coefficients obtained from regularized model at each \code{lambda}. This is motivated by formulating L0 variable selection in ADMM form, which shows significant improvement over the commonly used regularized methods without this technique.
}
\value{
An object with S3 class \code{"ADMMnet"}.
  \item{Beta}{a sparse Matrix of coefficients, stored in class "dgCMatrix".}
  \item{Beta0}{coefficients after additionally tuning the number of non-zeros, for \code{inzero = TRUE}.}
  \item{fit}{a data.frame containing \code{lambda} and the number of non-zero coefficients \code{nzero}. With cross-validation, additional results are reported, such as average cross-validation partial likelihood \code{cvm} and its standard error \code{cvse}, and \code{index} with `*' indicating the minimum \code{cvm}. For \code{family = "gaussian"}, \code{rsq} is also reported.}
  \item{fit0}{a data.frame containing \code{lambda}, \code{cvm} and \code{nzero} based on \code{inzero = TRUE}.}
  
  \item{lambda.min}{value of \code{lambda} that gives minimum \code{cvm}.}
  \item{lambda.opt}{value of \code{lambda} based on \code{inzero = TRUE}.}
  
  \item{penalty}{penalty type.}
  \item{adaptive}{logical flags for adaptive version (see above).}
  \item{flag}{convergence flag (for internal debugging). \code{flag = 0} means converged.}
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

\author{
Xiang Li, Shanghong Xie, Donglin Zeng and Yuanjia Wang\cr
Maintainer: Xiang Li <xl2473@cumc.columbia.edu>, Shanghong Xie <sx2168@cumc.columbia.edu>
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\section{Warning}{
  It may terminate and return \code{NULL}.
}

\seealso{
\code{\link{print.ADMMnet}}
}
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
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Number of non-zeros}
\keyword{Regularization}


