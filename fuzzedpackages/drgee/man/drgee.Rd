\encoding{latin1}
\name{drgee}
\alias{drgee}
\alias{coef.drgee}
\alias{vcov.drgee}
\alias{print.drgee}
\alias{summary.drgee}
\alias{print.summary.drgee}
%\alias{naiveVcov.drgee}
\alias{clusterRobustVcov.drgee}
\title{Doubly Robust Generalized Estimating Equations}
\description{
  \code{drgee} is used to estimate an exposure-outcome effect adjusted
  for  additional covariates. The estimation is based on regression
  models for the outcome, exposure or a combination of both.
  For clustered data the models may
  have cluster-specific intercepts.
}
\usage{
drgee(outcome, exposure,
      oformula, eformula, iaformula = formula(~1),
      olink = c("identity", "log", "logit"),
      elink = c("identity", "log", "logit"),
      data, subset = NULL, estimation.method = c("dr", "o", "e"),
      cond = FALSE, clusterid, clusterid.vcov, rootFinder = findRoots,
      intercept = TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{outcome}{
    The outcome as variable or as a character string naming a variable
    in the \code{data} argument. If missing, the outcome is assumed
    to be the response of \code{oformula}.
  }

  \item{exposure}{
    The exposure as variable or as a character string naming a variable
    in the \code{data} argument. If missing, the exposure is assumed
    to be the response of \code{eformula}.
  }

  \item{oformula}{
    An expression or formula for the outcome nuisance model.
  }

  \item{eformula}{
    An expression or formula for the exposure nuisance model.
  }

  \item{iaformula}{
    An expression or formula where the RHS should contain the variables
    that "interact" (i.e. are supposed to be multiplied with) with the
    exposure in the main model. "1" will always added. Default value is no
    interactions, i.e. \code{iaformula = formula(~1)}.
  }

  \item{olink}{
    A character string naming the link function in the outcome nuisance
    model. Has to be \code{"identity"}, \code{"log"} or
    \code{"logit"}.Default is \code{"identity"}.
  }

  \item{elink}{
    A character string naming the link function in the exposure nuisance
    model. Has to be \code{"identity"}, \code{"log"} or
    \code{"logit"}. Default is \code{"identity"}.
  }

  \item{data}{
    A data frame or environment containing the variables used.
    If missing, variables are expected to be found in the
    calling environment of the calling environment.
  }

  \item{subset}{
    An optional vector defining a subset of the data to be used.
  }

  \item{estimation.method}{
    A character string naming the desired estimation method. Choose
    \code{"o"} for O-estimation,
    \code{"e"} for E-estimation or
    \code{"dr"} for DR-estimation. Default is \code{"dr"}.
  }

  \item{cond}{
    A logical value indicating whether the nuisance models should have
    cluster-specific intercepts. Requires a \code{clusterid} argument.
  }

  \item{rootFinder}{
    A function to solve a system of non-linear equations. Default
    is \code{findRoots}.
  }

  \item{clusterid}{
    A cluster-defining variable or a character string naming a
    cluster-defining variable in the \code{data} argument. If it is not
    found in the \code{data} argument, it will be searched for in the
    calling frame. If missing, each observation will be considered to be
    a separate cluster. This argument is required when \code{cond = TRUE}.
  }

  \item{clusterid.vcov}{
    A cluster-defining variable or a character string naming a
    cluster-defining variable in the \code{data} argument to be used for
    adding contributions from the same cluster. These clusters can be
    different from the clusters defined by \code{clusterid}. However,
    each cluster defined by \code{clusterid} needs to be contained in
    exactly one cluster defined by \code{clusterid.vcov}. This variable
    is useful when the clusters are hierarchical.
  }

  \item{intercept}{
    A boolean to choose whether the nuisance parameters in doubly robust
    conditional logistic regression should be fitted with a model with an
    intercept. Only used for doubly robust condtional logistic regression.
  }
  
  \item{\dots}{
    Further arguments to be passed to the function \code{rootFinder}.
  }
}
\details{
  \code{drgee} estimates the parameter \eqn{\beta}{beta} in a main
  model \eqn{g\{E(Y|A,L)\}-g\{E(Y|A=0,L)\}=\beta^T \{A\cdot X(L)\}}{g\{E(Y|A,L)-g\{E(Y|A,L)\}=beta^T (A * X(L))\}},
  where \eqn{Y} is the outcome of interest, \eqn{A} is the exposure of
    interest, and \eqn{L} is a vector of covariates that we wish to
    adjust for. \eqn{X(L)} is a vector valued function of \eqn{L}. Note that \eqn{A
    \cdot X(L)}{A*X(L)} should be interpreted as a columnwise
  multiplication and that \eqn{X(L)} will always contain a column of 1's.
  Given a specification of an outcome nuisance model \eqn{g\{E(Y|A=0,L)=\gamma^T
    V(L)}{g\{E(Y|A=0,L)=gamma^T V(L)} (where \eqn{V(L)} is a function of \eqn{L})
  O-estimation is performed. Alternatively, leaving \eqn{g\{E(Y|A=0,L)}
  unspecified and using an exposure nuisance model \eqn{h\{E(A|L)\}=\alpha^T
    Z(L)}{h\{E(Y|L)\}=alpha^T Z(L)} (where \eqn{h} is a link
  function and \eqn{Z(L)} is a function of \eqn{L}) E-estimation
  is performed.  When \eqn{g} is logit, the exposure nuisance
  model is required be of the form
  \eqn{logit\{E(A|Y=0,L)\}=\alpha^T Z(L)}{logit\{E(A|Y=0,L)\}=alpha^T Z(L)}.
  In this case the exposure needs to binary.

  Given both an outcome and an exposure nuisance model, DR-estimation can be
  performed. DR-estimation gives a consistent estimate of the parameter
  \eqn{\beta}{beta} when either the outcome nuisance model or
  the exposure nuisance model
  is correctly specified, not necessarily both.

  Usage is best explained through an example. Suppose that we are
  interested in the parameter vector \eqn{(\beta_0,
    \beta_1)}{(beta_0,beta_1)} in a main model
  \eqn{logit\{E(Y|A,L_1,L_2)\}-logit\{E(Y|A=0,L_1,L_2)\}=\beta_0 A + \beta_1
    A \cdot L_1}{logit\{E(Y|A,L_1,L_2)-logit\{E(Y|A=0,L_1,L_2)\}=beta_0 A + beta_1 A *
    L_1} where \eqn{L_1} and \eqn{L_2} are the covariates that we wish
  to adjust for. To adjust for  \eqn{L_1} and \eqn{L_2}, we can use an outcome
  nuisance model \eqn{E(Y|A=0,L_1,L_2;\gamma_0, \gamma_1)=\gamma_0 + \gamma_1
  L_1}{E(Y|A=0,L_1,L_2;gamma_0, gamma_1)=gamma_0 + gamma_1 L_1} or an
  exposure nuisance model \eqn{logit\{E(A|Y=0,L_1,L_2)\}=\alpha_0+\alpha_1
  L_1+\alpha_2 L_2}{logit\{E(A|Y=0,L_1,L_2)\}=alpha_0+alpha_1
  L_1+alpha_2 L_2} to calculate estimates of \eqn{\beta_0}{beta_0} and \eqn{\beta_1}{beta_1}
  in the main model. We specify the outcome nuisance model as \code{oformula=Y~L_1}
  and \code{olink = "logit"}. The exposure nuisance model is specified as
  \code{eformula = A~L_1+L_2} and \code{elink = "logit"}.
  Since the outcome \eqn{Y} and the exposure \eqn{A} are
  identified as the LHS of \code{oformula} and \code{eformla}
  respectively and since the outcome link is specified in the
  \code{olink} argument,
  the only thing left to specify for the main model is the
  (multiplicative) interactions \eqn{A\cdot X(L)=A\cdot
  (1,L_1)^T}{A X(L)=A (1,L_1)^T}. This
  is done by specifying \eqn{X(L)} as
  \code{iaformula = ~L_1}, since \eqn{1} is always included in \eqn{X(L)}.
  We can then perform O-estimation, E-estimation or DR-estimation by
  setting \code{estimation.method} to \code{"o"},
  \code{"e"} or \code{"dr"} respectively. O-estimation uses only the
  outcome nuisance model, and E-estimation uses only the exposure
  nuisance model. DR-estimation uses both nuisance models, and gives a
  consistent estimate of \eqn{(\beta_0,\beta_1)}{(beta_0,beta_1)} if either nuisance model is correct, not necessarily both.

  When \code{estimation.method = "o"}, the RHS of \code{eformula} will be
  ignored. The \code{eformula} argument can also be replaced by an \code{exposure}
  argument specifying what the exposure of interest is.

  When \code{estimation.method = "e"}, the RHS of \code{oformula} will be
  ignored.  The \code{oformula} argument can also be replaced by an \code{outcome}
  argument specifying what the outcome of interest is.

  When \code{cond = TRUE} the nuisance models will be assumed to have
  cluster-specific intercept. These intercepts will not estimated.

  When E-estimation or DR-estimation is chosen with
  \code{olink = "logit"}, the exposure link will be
  changed to \code{"logit"}. Note that this choice
  of outcome link does not work for DR-estimation
  when \code{cond = TRUE}.

  Robust variance for the estimated parameter is calculated
  using the function \code{robVcov}. A cluster robust variance is calculated when
  a character string naming a cluster variable is
  supplied in the \code{clusterid} argument.

  For E-estimation when \code{cond = FALSE} and \eqn{g} is the identity
  or log link, see Robins et al. (1992).

  For DR-estimation when \code{cond = TRUE} and \eqn{g} is the identity
  or log link, see Robins (1999). For DR-estimation when
  \eqn{g} is the logit link, see Tchetgen et al. (2010).

  O-estimation can also be performed using the \code{gee} function.

}
\value{
  \code{drgee} returns an object of class \code{drgee} containing:

  \item{coefficients }{Estimates of the parameters in the main model.}
  \item{vcov }{Robust variance for all main model parameters.}
  \item{coefficients.all }{Estimates of all estimated parameters.}
  \item{vcov.all }{Robust variance of the all parameter estimates.}
  \item{optim.object }{An estimation object returned from the function specified
  	in the \code{rootFinder}, if this function is called for the
  estimation of the main model parameters.}
  \item{optim.object.o }{An estimation object returned from the function specified
  	in the \code{rootFinder}, if this function is called for the
  estimation of the outcome nuisance parameters.}
  \item{optim.object.e }{An estimation object returned from the function specified
  	in the \code{rootFinder}, if this function is called for the
  estimation of the outcome nuisance parameters.}
  \item{call }{The matched call.}
  \item{estimation.method }{The value of the input argument \code{estimation.method}.}
  \item{data}{The original data object, if given as an input argument}
  \item{oformula}{The original oformula object, if given as an input argument}
  \item{eformula}{The original eformula object, if given as an input argument}
  \item{iaformula}{The original iaformula object, if given as an input argument}

  The class methods \code{coef} and \code{vcov} can be used to extract
  the estimated parameters and their covariance matrix from a
  \code{drgee} object. \code{summary.drgee} produces a summary of the
  calculations.
}
\seealso{
  \code{\link{gee}} for O-estimation, \code{\link{findRoots}} for
  nonlinear equation solving and \code{\link{robVcov}} for
  estimation of variance.
}
\author{
Johan Zetterqvist, Arvid \enc{Sjölander}{Sjolander}
}

\examples{

## DR-estimation when
## the main model is
## E(Y|A,L1,L2)-E(Y|A=0,L1,L2)=beta0*A+beta1*A*L1
## and the outcome nuisance model is
## E(Y|A=0,L1,L2)=gamma0+gamma1*L1+gamma2*L2
## and the exposure nuisance model is
## E(A|Y=0,L1,L2)=expit(alpha0+alpha1*L1+alpha2*l2)

library(drgee)

expit<-function(x) exp(x)/(1+exp(x))

n<-5000

## nuisance
l1<-rnorm(n, mean = 0, sd = 1)
l2<-rnorm(n, mean = 0, sd = 1)

beta0<-1.5
beta1<-1
gamma0<--1
gamma1<--2
gamma2<-2
alpha0<-1
alpha1<-5
alpha2<-3

## Exposure generated from the exposure nuisance model
a<-rbinom(n,1,expit(alpha0 + alpha1*l1 + alpha2*l2))
## Outcome generated from the main model and the
## outcome nuisance model
y<-rnorm(n,
mean = beta0 * a + beta1 * a * l1 + gamma0 + gamma1 * l1 + gamma2 * l2,
sd = 1)

simdata<-data.frame(y,a,l1,l2)

## outcome nuisance model misspecified and
## exposure nuisance model correctly specified

## DR-estimation
dr.est <- drgee(oformula = formula(y~l1),
eformula = formula(a~l1+l2),
iaformula = formula(~l1),
olink = "identity", elink = "logit",
data = simdata, estimation.method = "dr")
summary(dr.est)

## O-estimation
o.est <- drgee(exposure = "a", oformula = formula(y~l1),
iaformula = formula(~l1), olink = "identity",
data = simdata, estimation.method = "o")
summary(o.est)

## E-estimation
e.est <- drgee(outcome = "y", eformula = formula(a~l1+l2),
iaformula = formula(~l1), elink="logit",
data = simdata, estimation.method = "e")
summary(e.est)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{models}
\keyword{regression}
\keyword{cluster}

\references{
  Orsini N., Belocco R., \enc{Sjölander}{Sjolander} A. (2013), Doubly
  Robust Estimation in Generalized Linear Models, \emph{Stata Journal},
  \bold{13}, 1, pp. 185--205

  Robins J.M., Mark S.D., Newey W.K. (1992), Estimating Exposure
  Effects by Modelling the Expectation of Exposure Conditional
  on Confounders, \emph{Biometrics}, \bold{48}, pp. 479--495

  Robins JM (1999), Robust Estimation in Sequentially Ignorable
  Missing Data and Causal Inference Models, \emph{Proceedings of the
    American Statistical Association Section on Bayesian Statistical
    Science}, pp. 6--10

  Tchetgen E.J.T., Robins J.M., Rotnitzky A. (2010), On Doubly Robust
  Estimation in a Semiparametric Odds Ratio Model, \emph{Biometrika},
  \bold{97}, 1, 171--180

  Zetterqvist J., Vansteelandt S., Pawitan Y.,
  \enc{Sjölander}{Sjolander} (2016), Doubly Robust Methods for Handling
  Confounding by Cluster, \emph{Biostatistics}, \bold{17}, 2, 264--276
}
