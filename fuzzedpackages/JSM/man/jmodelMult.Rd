
\name{jmodelMult}

\alias{jmodelMult}

\title{Semiparametric Joint Models for Survival and Longitudinal Data with Nonparametric Multiplicative Random Effects}

\description{
This function applies a maximum likelihood approach to fit the semiparametric joint models of survival and normal
longitudinal data. The survival model is assumed to come from a class of transformation models, including the Cox
proportional hazards model and the proportional odds model as special cases. The longitudinal process is modeled by  
nonparametric multiplicative random effects (NMRE) model.
}

\usage{
jmodelMult(fitLME, fitCOX, data, model = 1, rho = 0, timeVarY = NULL, 
           timeVarT = NULL, control = list(), \dots)
}

\arguments{
  \item{fitLME}{an object inheriting from class \code{lme} representing a fitted nonparametric multiplicative random effects model. See \bold{Details} and \bold{Note} and \bold{Examples}.}
  \item{fitCOX}{an object inheriting from class \code{coxph} representing a fitted Cox proportional hazards regression model. Specifying \code{x = TRUE} is required in the call to \code{coxph()} to include the design matrix in the object fit. See \bold{Note}.}
  \item{data}{a data.frame containing all the variables included in the joint modeling. See \bold{Note}.}
  \item{model}{an indicator specifying the dependency between the survival and longitudinal outcomes. Default is 1. See 
  \bold{Details}.}
  \item{rho}{a nonnegative real number specifying the transformation model you would like to fit. Default is 0, i.e. the Cox proportional hazards model. See \bold{Details}.}
  \item{timeVarY}{a character string indicating the time variable in the NMRE model. See \bold{Examples}.}
  \item{timeVarT}{a character string indicating the time variable in the \code{coxph} object. Normally it is \code{NULL}. See \bold{Note} and \bold{Examples}.}
  \item{control}{a list of control values for the estimation algorithm with components:
    \describe{
      \item{tol.P}{tolerance value for convergence in the parameters with default value 1e-03. See \bold{Details}.}
      \item{tol.L}{tolerance value for convergence in the log-likelihood with default value 1e-06. See \bold{Details}.}
      \item{max.iter}{the maximum number of EM iterations with default value 250.}
      \item{SE.method}{a character string specifying the standard error estimation method. Default is \code{"PRES"}. See \bold{Details} and \bold{Note}.}
      \item{delta}{a positive value used for numerical differentiation in the \code{SE.method}. Default is 1e-05 if \code{"PRES"} is used and 1e-03 otherwise. See \bold{Details}.}
      \item{nknot}{the number of Gauss-Hermite quadrature knots used to approximate the integrals over the random effects. Under the nonparametric multiplicative random effects model, there are only one-dimensional integrations and the default for \code{nknot} is 11.}
      }
  }
  \item{\dots}{additional options to be passed to the \code{control} argument.}
}
      
\details{
The \code{jmodelMult} function fits joint models for survival and longitudinal data. Nonparametric multiplicative random effects models (NMRE) are assumed for the longitudinal processes. With the Cox proportional hazards model and the proportional odds model as special cases, a general class of transformation models are assumed for the survival processes. The baseline hazard functions are left unspecified, i.e. no parametric forms are assumed, thus leading to semiparametric models. For detailed model formulation, please refer to Xu, Hadjipantelis and Wang (2017). 

The longitudinal model (NMRE) is written as \deqn{Y_i(t)=\mu_i(t) + \varepsilon_i(t) = b_i\times\mathbf{B}^\top(t)\boldsymbol\gamma + \varepsilon_i(t),} where \eqn{\mathbf{B}(t)=(B_1(t), \cdots, B_L(t))} is a vector of B-spline basis functions and \eqn{b_i} is a random effect \eqn{\sim\mathcal{N}(1, \sigma_b^2)}. Note that we also allow the inclusion of baseline covariates as columns of \eqn{\mathbf{B}(t)}. If \code{model = 1}, then the linear predictor for the survival model is expressed as \deqn{\eta(t) = \mathbf{W}_i^\top(t)\boldsymbol\phi + \alpha\mu_i(t),} indicating that the entire longitudinal process (free of error) enters the survival model as a covariate. If other values are assigned to the \code{model} argument, the linear predictor for the surival model is then expressed as \deqn{\eta(t) = \mathbf{W}_i^\top(t)\boldsymbol\phi + \alpha b_i,} suggesting that the survival and longitudinal models only share the same random effect.

The survival model is written as \deqn{\Lambda(t|\eta(t)) = G\left[\int_0^t\exp\{\eta(s)\}d\Lambda_0(s)\right],} where \eqn{G(x) = \log(1 + \rho x) / \rho} with \eqn{\rho \geq 0} is the class of logarithmic transfomrations. If \code{rho = 0}, then \eqn{G(x) = x}, yielding the Cox proportional hazards model. If \code{rho = 1}, then \eqn{G(x) = \log(1 + x)}, yielding the proportional odds model. Users could assign any nonnegative real value to \code{rho}.

An expectation-maximization (EM) algorithm is implemented to obtain parameter estimates. The convergence criterion is either of (i) \eqn{\max \{ | \boldsymbol\theta^{(t)} - \boldsymbol\theta^{(t - 1)} | / ( | \boldsymbol\theta^{(t - 1)} | + .Machine\$double.eps \times 2 ) \} < tol.P}, or (ii) \eqn{| L(\boldsymbol\theta^{(t)}) - L(\boldsymbol\theta^{(t - 1)})| / ( | L(\theta^{(t - 1)}) | + .Machine\$double.eps \times 2 ) < tol.L}, is satisfied. Here \eqn{\boldsymbol\theta^{(t)}} and \eqn{\boldsymbol\theta^{(t-1)}} are the vector of parameter estimates at the \eqn{t}-th and \eqn{(t-1)}-th EM iterations, respectively; \eqn{L(\boldsymbol\theta)} is the value of the log-likelihood function evaluated at \eqn{\boldsymbol\theta}. Users could specify the tolerance values \code{tol.P} and \code{tol.L} through the \code{control} argument.

For standard error estimation for the parameter estimates, three methods are provided, namely \code{"PRES"}, \code{"PFDS"} and \code{"PLFD"} (detailed information are referred to Xu, Baines and Wang (2014)). In the \code{control} argument, if \code{SE.method = "PRES"}, numerically differentiating the profile Fisher score vector with Richardson extrapolation is applied; if \code{SE.method = "PFDS"}, numerically differentiating the profile Fisher score vector with forward difference is applied; if \code{SE.method = "PLFD"}, numerially (second) differentiating the profile likelihood with forward difference is applied. Generally, numerically differentiating a function \eqn{f(x)} (an arbitrary function) with forward difference is expressed as \deqn{f^\prime(x) = \frac{f(x + \delta) - f(x)}{\delta},} and that with Richardson extrapolation is expressed as \deqn{f^\prime(x) = \frac{f(x - 2\delta) - 8f(x - \delta) + 8f(x + \delta) - f(x + 2\delta)}{12\delta}.} Users could specify the value of \eqn{\delta} through the \code{delta} item in the \code{control} argument.
}

\value{
See \code{\link{jmodelMultObject}} for the components of the fit.
}

\references{
Dabrowska, D. M. and Doksun K. A. (1988) Partial Likelihood in Transformation Models with Censored Data. \emph{Scandinavian Journal of Statistics} \bold{15}, 1--23.

Ding, J. and Wang, J. L. (2008) Modeling longitudinal data with nonparametric multiplicative random effects jointly with survival data. \emph{Biometrics} \bold{64}, 546--556.

Tsiatis, A. A. and Davidian, M. (2004) Joint modeling of longitudinal and time-to-event data: an overview. 
\emph{Statistica Sinica} \bold{14}, 809--834.

Xu, C., Baines, P. D. and Wang, J. L. (2014) Standard error estimation using the EM algorithm for the joint modeling of survival and longitudinal data. \emph{Biostatistics} \bold{15}, 731--744

Xu, C., Hadjipantelis, P. Z. and Wang, J. L. (2020) Semiparametric joint modeling of survival and longitudinal data: the R package JSM. \emph{Journal of Statistical Software} <doi:10.18637/jss.v093.i02>.

Zeng, D. and Lin, D. (2007) Maximum likelihood estimation in semiparametric regression models with censored data. \emph{Journal of the Royal Statistical Society: Series B} \bold{69}, 507--564.
}

\author{
Cong Xu \email{helenxu1112@gmail.com}
Pantelis Z. Hadjipantelis \email{pantelis@ucdavis.edu}
}

\note{
1. To fit a nonparametric multiplicative random effects model, the fixed effect in the \code{fitLME} object should be a matrix of B-spline basis functions (an object from the \code{bs} function) with the possibility of including baseline covariates and the random effect should only include a random intercept. In the \code{bs} function, it is a good practice to specify the boundary knots through the \code{Boundary.knots} argument, where the upper boundary knot is typically the longest follow-up time among all subjects. See \bold{Examples}.

2. Currently, \code{jmodelMult()} could only handle the \code{fitLME} object with a simple random-effects structure (only the \code{pdDiag()} class). Moreover, the within-group correlation and heteroscedasticity structures in the \code{fitLME} object (i.e. the \code{correlation} and \code{weights} argument of \code{lme()}) are ignored.

3. The \code{data} argument in \code{jmodelMult()}, \code{lme()} and \code{coxph()} should be the same data frame.

4. For the \code{fitCOX} object, only the \eqn{W_i(t)} in the linear predictor \eqn{\eta(t)} for the survial model (see \bold{Details}) should be involved in the \code{formula} argument of \code{coxph{}}. Since \code{coxph()} uses the same data frame as \code{lme()} does, a time-dependent Cox model must be fitted by \code{coxph()} although \eqn{W_i(t)} may only contain time-independent covariates. See \bold{Examples}.

5. If \eqn{W_i(t)} in the linear predictor \eqn{\eta(t)} for the survial model (see \bold{Details}) does involve time-dependent covariate, then \code{timeVarT} must specify the name of the time variable involved (see \bold{Examples}).

6. The standard error estimates are obtained by numerical approximations which is naturally subject to numerical errors. Therefore, in extreme cases, there may be \code{NA} values for part of the standard error estimates.
}

\seealso{
\code{\link{jmodelMultObject}},
\code{\link{lme}},
\code{\link{coxph}},
\code{\link{Surv}},
\code{\link{bs}}
}

\examples{
# linear mixed-effects model fit where the fixed effect is modeled by 
# quadratic B-splie basis with no internal knots
fitLME <- lme(log(serBilir) ~ bs(obstime, degree = 2, Boundary.knots = c(0, 15)), 
              random = ~ 1 | ID, data = pbc)
# Cox proportional hazards model fit with a single time-independent covariate
fitCOX <- coxph(Surv(start, stop, event) ~ drug, data = pbc, x = TRUE)

# joint model fit which assumes the Cox proportional hazards model for the survival process 
# and NMRE for the longitudinal process. Use 'max.iter = 25', 'nknot = 3' and 
# the 'PFDS' method to calculate standard error estimates as a quick toy example 
fitJTMult.ph <- jmodelMult(fitLME, fitCOX, pbc, timeVarY = "obstime", 
                           control = list(SE.method = 'PFDS', max.iter = 25, nknot = 3))
summary(fitJTMult.ph)

\dontrun{
# joint model fit with the default control
fitJTMult.ph2 <- jmodelMult(fitLME, fitCOX, pbc, timeVarY = "obstime")
summary(fitJTMult.ph2)
# joint model fit where the survival and longitudinal processes only share 
# the same random effect
fitJTMult.ph3 <- jmodelMult(fitLME, fitCOX, pbc, model = 2, timeVarY = "obstime")
summary(fitJTMult.ph3)

# joint model fit which assumes the proportional odds model for the survival process 
# and NMRE model for the longitudinal process
fitJTMult.po <- jmodelMult(fitLME, fitCOX, pbc, rho = 1, timeVarY = "obstime")
summary(fitJTMult.po)
# joint model fit where the survival and longitudinal processes only share 
# the same random effect
fitJTMult.po2 <- jmodelMult(fitLME, fitCOX, pbc, model = 2, rho = 1, timeVarY = "obstime")
summary(fitJTMult.po2)

# allow baseline covariates in the NMRE model for the longitudinal process
fitLME2 <- lme(log(serBilir) ~ drug + bs(obstime, degree = 2, Boundary.knots = c(0, 15)), 
               random = ~1 | ID, data = pbc)
fitJTMult.ph4 <- jmodelMult(fitLME2, fitCOX, pbc, timeVarY = "obstime")
summary(fitJTMult.ph4)

# Cox proportional hazards model fit with a time-dependent covariate
fitCOX2 <- coxph(Surv(start, stop, event) ~ drug + as.numeric(drug) : obstime, 
                 data = pbc, x = TRUE)
# joint model fit in which \code{timeVarT} must be specified
fitJTMult.ph5 <- jmodelMult(fitLME, fitCOX2, pbc, timeVarY = "obstime", timeVarT = 'obstime', 
                            control = list(max.iter = 300))
summary(fitJTMult.ph5)
}
}

\keyword{survival}
\keyword{models}
