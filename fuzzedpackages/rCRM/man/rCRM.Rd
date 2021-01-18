\name{rCRM}
\alias{rCRM}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Fit a 2-parameter CRM
}
\description{
Fit a 2-parameter CRM model (O'Quigley and Shen 1996) regularized with L2 norm (Friedman et al. 2010) adjusted by the distance with the target DLT rate.

}
\usage{
rCRM(x, y, dose0, tp = 0.3, family = "2P", mlambda = 1, nlambda = 50, rlambda = NULL,
wldose = NULL, nfolds = length(y), foldid = NULL, keep.beta = FALSE,
thresh = 1e-07, maxit = 1e+04, threshP = 1e-06, threshB = 100)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{input vector of dose.
}
  \item{y}{response variable. \code{y} is a binary vector with 0 (not DLT) and 1 (DLT).
}
  \item{dose0}{dose regimen. \code{x} should be included in \code{dose0}.
}
  \item{tp}{target toxicity probability. Default is 0.3.
}
  \item{family}{type of CRM models. Now only supports 2-paraemter CRM, \code{2P}.
}
  \item{mlambda}{maximum of tuning parameter \code{lambda}. The optimal \code{lambda} is selected by cross-validation.
}
  \item{nlambda}{number of \code{lambda} values. Default is 50.
}
  \item{rlambda}{fraction of \code{mlambda} to determine the smallest value for \code{lambda}. The default is \code{rlambda = 0.0001} when the number of observations is larger than or equal to the number of variables; otherwise, \code{rlambda = 0.01}.
}
  \item{wldose}{penalty weights used with L2 norm (adaptive L2). The \code{wldose} is a vector of non-negative values with the same length as \code{dose0}. Default is NULL indicating that weights are calculated based on MLE.
}
  \item{nfolds}{number of folds. With \code{nfolds = 1} and \code{foldid = NULL}, cross-validation is not performed. For cross-validation, smallest value allowable is \code{nfolds = 3}. Specifying \code{foldid} overrides \code{nfolds}. Default is \code{nfolds=length(y)} indicating leave-one-out cross-validation.
}
  \item{foldid}{an optional vector of values between 1 and \code{nfolds} specifying which fold each observation is in. Default is \code{foldid=NULL}.
}
  \item{keep.beta}{logical flag for returning estimates for all \code{lambda} values. For \code{keep.beta = FALSE}, only return the estimate with the minimum cross-validation value.
}
  \item{thresh}{convergence threshold for coordinate descent. Default value is \code{1E-7}.
}
  \item{maxit}{maximum number of iterations for coordinate descent. Default is \code{1E+4}.
}
  \item{threshP}{boundary for calculating the probability of DLT. Default is \code{1E-6}. The estimated probability is truncated between \code{1E-6} and 1-\code{1E-6}.
}
  \item{threshB}{boundary for calculating the parameters. Default is 100.  The estimates are truncated between \code{-100} and \code{100}.
}
}
\details{
 One-step coordinate descent algorithm is applied for each \code{lambda}.
 Cross-validation is used for tuning parameters.
}
\value{
An object with S3 class \code{"rCRM"}.
  \item{Beta}{estimates in 2-parameter CRM model.}
  \item{fit}{a data.frame containing \code{lambda} and proportion of deviance. With cross-validation, additional results are reported, such as average cross-validation likelihood \code{cvm} and its standard error \code{cvse}, and \code{index} with `*' indicating the minimum \code{cvm}.}


  \item{lambda.min}{value of \code{lambda} that gives minimum \code{cvm}.}
  \item{flag}{convergence flag (for internal debugging). \code{flag = 0} means converged.}
  \item{prob}{estimated probability of DLT at each \code{dose0}.}
  \item{dose.close}{the index of dose in \code{dose0} with the \code{prob} cloest to \code{tp}.}
  \item{family}{type of CRM models. \code{2P} is 2-parameter CRM model.}
}

\references{O'Quigley, J., Shen, L.Z. (1996).
\emph{Continual reassessment method: a likelihood approach. Biometrics, 673-684.}\cr
Friedman, J., Hastie, T. and Tibshirani, R. (2010).
\emph{Regularization paths for generalized linear models via coordinate descent, Journal of Statistical Software, Vol. 33(1), 1.}\cr
}

\author{
Maintainer: Xiang Li <xli256@its.jnj.com>
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\section{Warning}{
  It may terminate and return \code{NULL}.
}


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
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{CRM}% use one of  RShowDoc("KEYWORDS")
\keyword{Regularization}% __ONLY ONE__ keyword per line

