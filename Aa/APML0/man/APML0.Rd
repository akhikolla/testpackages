\name{APML0}
\alias{APML0}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Fit a Model with Various Regularization Forms
}
\description{
Fit linear, logistic and Cox models regularized with L0, lasso (L1), elastic-net (L1 and L2), or net (L1 and Laplacian) penalty, and their adaptive forms, such as adaptive lasso / elastic-net and net adjusting for signs of linked coefficients.

It solves L0 penalty problem by simultaneously selecting regularization parameters and performing hard-thresholding (or selecting number of non-zeros). This augmented and penalized minimization method provides an approximation solution to the L0 penalty problem and runs as fast as L1 regularization.

The function uses one-step coordinate descent algorithm and runs extremely fast by taking into account the sparsity structure of coefficients. It could deal with very high dimensional data.
}
\usage{
APML0(x, y, family=c("gaussian", "binomial", "cox"), penalty=c("Lasso","Enet", "Net"),
Omega=NULL, alpha=1.0, lambda=NULL, nlambda=50, rlambda=NULL, wbeta=rep(1,ncol(x)),
sgn=rep(1,ncol(x)), nfolds=1, foldid=NULL, ill=TRUE, iL0=TRUE, icutB=FALSE, ncutB=10,
ifast=TRUE, isd=FALSE, iysd=FALSE, ifastr=TRUE, keep.beta=FALSE,
thresh=1e-6, maxit=1e+5, threshC=1e-5, maxitC=1e+2, threshP=1e-5)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{input matrix. Each row is an observation vector.
}
  \item{y}{response variable. For \code{family = "gaussian"},
  \code{y} is a continuous vector. For \code{family = "binomial"},
  \code{y} is a binary vector with 0 and 1. For \code{family = "cox"}, \code{y} is a two-column matrix with columns named `time' and `status'. `status' is a binary variable, with `1' indicating event, and `0' indicating right censored.
}
  \item{family}{type of outcome. Can be "gaussian", "binomial" or "cox".
}
  \item{penalty}{penalty type. Can choose \code{"Net"}, \code{"Enet"} (elastic net) and \code{"Lasso"}. For \code{"Net"}, need to specify \code{Omega}; otherwise, \code{"Enet"} is performed. For \code{penalty = "Net"}, the penalty is defined as \deqn{\lambda*{\alpha*||\beta||_1+(1-\alpha)/2*(\beta^{T}L\beta)},}
where \eqn{L} is a Laplacian matrix calculated from \code{Omega}.
}
  \item{Omega}{adjacency matrix with zero diagonal and non-negative off-diagonal, used for \code{penalty = "Net"} to calculate Laplacian matrix.
}
  \item{alpha}{ratio between L1 and Laplacian for \code{"Net"}, or between L1 and L2 for \code{"Enet"}. Default is \code{alpha = 1.0}, i.e. lasso.
}
  \item{lambda}{a user supplied decreasing sequence. If \code{lambda = NULL}, a sequence of \code{lambda} is generated based on \code{nlambda} and \code{rlambda}. Supplying a value of \code{lambda} overrides this.
}
  \item{nlambda}{number of \code{lambda} values. Default is 50.
}
  \item{rlambda}{fraction of \code{lambda.max} to determine the smallest value for \code{lambda}. The default is \code{rlambda = 0.0001} when the number of observations is larger than or equal to the number of variables; otherwise, \code{rlambda = 0.01}.
}
  \item{wbeta}{penalty weights used with L1 penalty (adaptive L1), given by \eqn{\sum_{j=1}^qw_j|\beta_j|.} The \code{wbeta} is a vector of non-negative values and works as adaptive L1. No penalty is imposed for those coefficients with zero values in \code{wbeta}. Default is 1 for all coefficients. The same weights are also applied to L0.
}
  \item{sgn}{sign adjustment used with Laplacian penalty (adaptive Laplacian). The \code{sgn} is a vector of 1 or -1. The \code{sgn} could be based on an initial estimate of \eqn{\beta}, and 1 is used for \eqn{\beta>0} and -1 is for \eqn{\beta<0}. Default is 1 for all coefficients.
}
  \item{nfolds}{number of folds. With \code{nfolds = 1} and \code{foldid = NULL} by default, cross-validation is not performed. For cross-validation, smallest value allowable is \code{nfolds = 3}. Specifying \code{foldid} overrides \code{nfolds}.
}
  \item{foldid}{an optional vector of values between 1 and \code{nfolds} specifying which fold each observation is in.
}
  \item{ill}{logical flag for using likelihood-based as the cross-validation criteria. Default is \code{ill = TRUE}. For \code{family = "gaussian"}, set \code{ill = FALSE} to use predict mean squared error as the criteria.
}
  \item{iL0}{logical flag for simultaneously performing L0-norm via performing hard-thresholding or selecting number of non-zeros. Default is \code{iL0 = TRUE}.
}
  \item{icutB}{logical flag for performing hard-thresholding by selecting the number of non-zero coefficients with the default of \code{icutB = FALSE}. Alternative way is to apply thresholding on the coefficients by setting \code{icutB = TRUE}.
}
  \item{ncutB}{the number of thresholds used for \code{icutB = TRUE}. Default is \code{ncutB=10}. Increasing \code{ncutB} may improve the variable selection performance but will increase the computation time.
}
  \item{ifast}{logical flag for searching for the best cutoff or the number of non-zero. Default is \code{ifast=TRUE} for local searching. Setting \code{ifast=TRUE} will search from the smallest cutoff (or number of non-zeros) to the largest but will increase the computation time.
}

  \item{isd}{logical flag for outputting standardized coefficients. \code{x} is always standardized prior to fitting the model. Default is \code{isd = FALSE}, returning \eqn{\beta} on the original scale.
}
  \item{iysd}{logical flag for standardizing \code{y} prior to computation, for \code{family = "gaussian"}. The returning coefficients are always based the original \code{y} (unstandardized). Default is \code{isd = FALSE}.
}
  \item{ifastr}{logical flag for efficient calculation of risk set updates for \code{family = "cox"}. Default is \code{ifastr = TRUE}. Setting \code{ifastr = FALSE} may improve the accuracy of calculating the risk set.
}
  \item{keep.beta}{logical flag for returning estimates for all \code{lambda} values. For \code{keep.beta = FALSE}, only return the estimate with the minimum cross-validation value.
}
  \item{thresh}{convergence threshold for coordinate descent. Default value is \code{1E-6}.
}
  \item{maxit}{Maximum number of iterations for coordinate descent. Default is \code{10^5}.
}
  \item{threshC}{convergence threshold for hard-thresholding for \code{family = "binomial"}. Default value is \code{1E-5}.
}
  \item{maxitC}{Maximum number of iterations for hard-thresholding for \code{family = "binomial"}. Default is \code{10^2}.
}
  \item{threshP}{Cutoff when calculating the probability in \code{family = "binomial"}. The probability is bounded within \code{threshP} and \code{1-threshP}. Default value is \code{1E-5}.
}
}
\details{
  One-step coordinate descent algorithm is applied for each \code{lambda}. Cross-validation is used for tuning parameters. For \code{iL0 =  TRUE}, we further perform hard-thresholding (for \code{icutB=TRUE}) to the coefficients or select the number of non-zero coefficients (for \code{icutB=FALSE}), which is obtained from regularized model at each \code{lambda}. This is motivated by formulating L0 variable selection in an augmented form, which shows significant improvement over the commonly used regularized methods without this technique. Details could be found in our publication.

  \code{x} is always standardized prior to fitting the model and the estimate is returned on the original scale for \code{isd=FALSE}.

  Each one element of \code{wbeta} corresponds to each variable in \code{x}. Setting the value in \code{wbeta} will not impose any penalty on that variable.

  For \code{family = "cox"}, \code{ifastr = TRUE} adopts an efficient way to update risk set and sometimes the algorithm ends before all \code{nlambda} values of \code{lambda} have been evaluated. To evaluate small values of \code{lambda}, use \code{ifast = FALSE}. The two methods only affect the efficiency of algorithm, not the estimates.

  \code{ifast = TRUE} seems to perform well.
}
\value{
An object with S3 class \code{"APML0"}.

  \item{a}{the intercept for \code{family = "gaussian"}.}
  \item{Beta}{a sparse Matrix of coefficients, stored in class "dgCMatrix". For \code{family = "binomial"}, the first coefficient is the intercept.}
  \item{Beta0}{coefficients after additionally performing L0-norm for \code{iL0 = TRUE}. For \code{family = "binomial"}, the first coefficient is the intercept.}
  \item{fit}{a data.frame containing \code{lambda} and the number of non-zero coefficients \code{nzero}. With cross-validation, additional results are reported, such as average cross-validation partial likelihood \code{cvm} and its standard error \code{cvse}, and \code{index} with `*' indicating the minimum \code{cvm}. For \code{family = "gaussian"}, \code{rsq} is also reported.}
  \item{fit0}{a data.frame containing \code{lambda}, \code{cvm} and \code{nzero} based on \code{iL0 = TRUE}. \code{cvm} in \code{fit0} may be different from \code{cvm} in \code{fit}, because the constaint on the number of non-zeros is imposed in the cross-validation. The maximum number of non-zeros is based on the full dataset not the one used in the cross-validation.}

  \item{lambda.min}{value of \code{lambda} that gives minimum \code{cvm}.}
  \item{lambda.opt}{value of \code{lambda} based on \code{iL0 = TRUE}.}

  \item{penalty}{penalty type.}
  \item{adaptive}{logical flags for adaptive version (see above).}
  \item{flag}{convergence flag (for internal debugging). \code{flag = 0} means converged.}
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

\author{
Xiang Li, Shanghong Xie, Donglin Zeng and Yuanjia Wang\cr
Maintainer: Xiang Li <spiritcoke@gmail.com>
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\section{Warning}{
  It may terminate and return \code{NULL}.
}

\seealso{
\code{\link{APML0}}, \code{\link{print.APML0}}
}
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
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{L0}
\keyword{Hard-thresholding}
\keyword{Number of non-zero}
\keyword{Regularization}


