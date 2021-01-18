
\name{jmodelTM}

\alias{jmodelTM}

\title{Semiparametric Joint Models for Survival and Longitudinal Data}

\description{
This function applies a maximum likelihood approach to fit the semiparametric joint models of survival and normal
longitudinal data. The survival model is assumed to come from a class of transformation models, including the Cox
proportional hazards model and the proportional odds model as special cases. The longitudinal process is modeled by  
liner mixed-effects models.
}

\usage{
jmodelTM(fitLME, fitCOX, data, model = 1, rho = 0, timeVarY = NULL, 
         timeVarT = NULL, control = list(), \dots)
}

\arguments{
  \item{fitLME}{an object inheriting from class \code{lme} representing a fitted linear mixed-effects model. See \bold{Note}.}
  \item{fitCOX}{an object inheriting from class \code{coxph} representing a fitted Cox proportional hazards regression model. Specifying \code{x = TRUE} is required in the call to \code{coxph()} to include the design matrix in the object fit. See \bold{Note}.}
  \item{data}{a data.frame containing all the variables included in the joint modeling. See \bold{Note}.}
  \item{model}{an indicator specifying the dependency between the survival and longitudinal outcomes. Default is 1. See 
  \bold{Details}.}
  \item{rho}{a nonnegative real number specifying the transformation model you would like to fit. Default is 0, i.e. the Cox proportional hazards model. See \bold{Details}.}
  \item{timeVarY}{a character string indicating the time variable in the linear mixed-effects model. See \bold{Examples}.}
  \item{timeVarT}{a character string indicating the time variable in the \code{coxph} object. Normally it is \code{NULL}. See \bold{Note} and \bold{Examples}.}
  \item{control}{a list of control values for the estimation algorithm with components:
    \describe{
      \item{tol.P}{tolerance value for convergence in the parameters with default value 1e-03. See \bold{Details}.}
      \item{tol.L}{tolerance value for convergence in the log-likelihood with default value 1e-06. See \bold{Details}.}
      \item{max.iter}{the maximum number of EM iterations with default value 250.}
      \item{SE.method}{a character string specifying the standard error estimation method. Default is \code{"PRES"}. See \bold{Details} and \bold{Note}.}
      \item{delta}{a positive value used for numerical differentiation in the \code{SE.method}. Default is 1e-05 if \code{"PRES"} is used and 1e-03 otherwise. See \bold{Details}.}
      \item{nknot}{the number of Gauss-Hermite quadrature knots used to approximate the integrals over the random effects. Default is 9 and 7 for one- and two-dimensional integration, respectively, and 5 for those with higher dimensions.}
    }
  }
  \item{\dots}{additional options to be passed to the \code{control} argument.}
}

\details{
The \code{jmodelTM} function fits joint models for survival and longitudinal data. Linear mixed-effects models are assumed for the longitudinal processes. With the Cox proportional hazards model and the proportional odds model as special cases, a general class of transformation models are assumed for the survival processes. The baseline hazard functions are left unspecified, i.e. no parametric forms are assumed, thus leading to semiparametric models. For detailed model formulation, please refer to Xu, Baines and Wang (2014). 

The longitudinal model is written as \deqn{Y_i(t)=\mu_i(t) + \varepsilon_i(t) = \mathbf{X}_i^\top(t)\boldsymbol\beta + \mathbf{Z}_i^\top(t)\mathbf{b}_i + \varepsilon_i(t).}, then the linear predictor for the survival model is expressed as \deqn{\eta(t) = \mathbf{W}_i^\top(t)\boldsymbol\phi + \alpha\mu_i(t),} indicating that the entire longitudinal process (free of error) enters the survival model as a covariate. If other values are assigned to the \code{model} argument, the linear predictor for the surival model is then expressed as \deqn{\eta(t) = \mathbf{W}_i^\top(t)\boldsymbol\phi + \alpha\mathbf{Z}_i^\top(t)\mathbf{b}_i,} suggesting that the survival and longitudinal models only share the same random effects.

The survival model is written as \deqn{\Lambda(t|\eta(t)) = G\left[\int_0^t\exp{\eta(s)}d\Lambda_0(s)\right],} where \eqn{G(x) = \log(1 + \rho x) / \rho} with \eqn{\rho \geq 0} is the class of logarithmic transfomrations. If \code{rho = 0}, then \eqn{G(x) = x}, yielding the Cox proportional hazards model. If \code{rho = 1}, then \eqn{G(x) = \log(1 + x)}, yielding the proportional odds model. Users could assign any nonnegative real value to \code{rho}. 

An expectation-maximization (EM) algorithm is implemented to obtain parameter estimates. The convergence criterion is either of (i) \eqn{\max \{ | \boldsymbol\theta^{(t)} - \boldsymbol\theta^{(t - 1)} | / ( | \boldsymbol\theta^{(t - 1)} | + .Machine\$double.eps \times 2 ) \} < tol.P}, or (ii) \eqn{| L(\boldsymbol\theta^{(t)}) - L(\boldsymbol\theta^{(t - 1)})| / ( | L(\theta^{(t - 1)}) | + .Machine\$double.eps \times 2 ) < tol.L}, is satisfied. Here \eqn{\boldsymbol\theta^{(t)}} and \eqn{\boldsymbol\theta^{(t-1)}} are the vector of parameter estimates at the \eqn{t}-th and \eqn{(t-1)}-th EM iterations, respectively; \eqn{L(\boldsymbol\theta)} is the value of the log-likelihood function evaluated at \eqn{\boldsymbol\theta}. Users could specify the tolerance values \code{tol.P} and \code{tol.L} through the \code{control} argument.

For standard error estimation for the parameter estimates, three methods are provided, namely \code{"PRES"}, \code{"PFDS"} and \code{"PLFD"} (detailed information are referred to Xu, Baines and Wang (2014)). In the \code{control} argument, if \code{SE.method = "PRES"}, numerically differentiating the profile Fisher score vector with Richardson extrapolation is applied; if \code{SE.method = "PFDS"}, numerically differentiating the profile Fisher score vector with forward difference is applied; if \code{SE.method = "PLFD"}, numerially (second) differentiating the profile likelihood with forward difference is applied. Generally, numerically differentiating a function \eqn{f(x)} (an arbitrary function) with forward difference is expressed as \deqn{f^\prime(x) = \frac{f(x + \delta) - f(x)}{\delta},} and that with Richardson extrapolation is expressed as \deqn{f^\prime(x) = \frac{f(x - 2\delta) - 8f(x - \delta) + 8f(x + \delta) - f(x + 2\delta)}{12\delta}.} Users could specify the value of \eqn{\delta} through the \code{delta} item in the \code{control} argument.
}

\value{
See \code{\link{jmodelTMObject}} for the components of the fit.
}

\references{
Dabrowska, D. M. and Doksun K. A. (1988) Partial Likelihood in Transformation Models with Censored Data. \emph{Scandinavian Journal of Statistics} \bold{15}, 1--23.

Tsiatis, A. A. and Davidian, M. (2004) Joint modeling of longitudinal and time-to-event data: an overview. \emph{Statistica Sinica} \bold{14}, 809--834.

Wulfsohn, M. S. and Tsiatis, A. A. (1997) A joint model for survival and longitudinal data measured with error. \emph{Biometrics} \bold{53}, 330--339.

Xu, C., Baines, P. D. and Wang, J. L. (2014) Standard error estimation using the EM algorithm for the joint modeling of survival and longitudinal data. \emph{Biostatistics} \bold{15}, 731--744.

Xu, C., Hadjipantelis, P. Z. and Wang, J. L. (2020) Semiparametric joint modeling of survival and longitudinal data: the R package JSM. \emph{Journal of Statistical Software} <doi:10.18637/jss.v093.i02>.

Zeng, D. and Lin, D. (2007) Maximum likelihood estimation in semiparametric regression models with censored data. \emph{Journal of the Royal Statistical Society: Series B} \bold{69}, 507--564.
}

\author{
Cong Xu \email{helenxu1112@gmail.com}
Pantelis Z. Hadjipantelis \email{pantelis@ucdavis.edu}
}

\note{
1. Currently, \code{jmodelTM()} could only handle the \code{fitLME} object with a simple random-effects structure (only the \code{pdDiag()} class). Moreover, the within-group correlation and heteroscedasticity structures in the \code{fitLME} object (i.e. the \code{correlation} and \code{weights} argument of \code{lme()}) are ignored.

2. The \code{data} argument in \code{jmodelTM()}, \code{lme()} and \code{coxph()} should be the same data frame.

3. For the \code{fitCOX} object, only the \eqn{W_i(t)} in the linear predictor \eqn{\eta(t)} for the survial model (see \bold{Details}) should be involved in the \code{formula} argument of \code{coxph{}}. Since \code{coxph()} uses the same data frame as \code{lme()} does, a time-dependent Cox model must be fitted by \code{coxph()} although \eqn{W_i(t)} may only contain time-independent covariates. See \bold{Examples}.

4. If \eqn{W_i(t)} in the linear predictor \eqn{\eta(t)} for the survial model (see \bold{Details}) does involve time-dependent covariate, then \code{timeVarT} must specify the name of the time variable involved. See \bold{Examples}.

5. The standard error estimates are obtained by numerical approximations which is naturally subject to numerical errors. Therefore, in extreme cases, there may be \code{NA} values for part of the standard error estimates.
}

\seealso{
\code{\link{jmodelTMObject}},
\code{\link{lme}},
\code{\link{coxph}},
\code{\link{Surv}}
}

\examples{
# linear mixed-effects model fit with random intercept
fitLME <- lme(sqrt(CD4) ~ obstime + I(obstime ^ 2) + drug : obstime + drug : I(obstime ^ 2), 
              random = ~ 1 | ID, data = aids)
# Cox proportional hazards model fit with a single time-independent covariate
fitCOX <- coxph(Surv(start, stop, event) ~ drug, data = aids, x = TRUE)

# joint model fit which assumes the Cox proportional hazards model for the survival process
# Use 'max.iter = 5', 'nknot = 3' and the 'PFDS' method to calculate standard 
# error estimates as a quick toy example 
fitJT.ph <- jmodelTM(fitLME, fitCOX, aids, timeVarY = 'obstime', 
                     control = list(SE.method = 'PFDS', max.iter = 5, nknot = 3))
summary(fitJT.ph)

\dontrun{
# joint model fit with the default control
fitJT.ph2 <- jmodelTM(fitLME, fitCOX, aids, timeVarY = 'obstime')
summary(fitJT.ph2)
# joint model fit where the survival and longitudinal processes only share 
# the same random effect
fitJT.ph3 <- jmodelTM(fitLME, fitCOX, aids, model = 2, timeVarY = 'obstime')
summary(fitJT.ph3)

# joint model fit which assumes the proportional odds model for the survival process
fitJT.po <- jmodelTM(fitLME, fitCOX, aids, rho = 1, timeVarY = 'obstime')
summary(fitJT.po)
# joint model fit where the survival and longitudinal processes only share 
# the same random effect
fitJT.po2 <- jmodelTM(fitLME, fitCOX, aids, model = 2, rho = 1, timeVarY = 'obstime')
summary(fitJT.po2)

# linear mixed-effects model fit with random intercept and random slope
fitLME2 <- lme(sqrt(CD4) ~ drug + obstime + I(obstime ^ 2) + drug : obstime + 
               drug : I(obstime ^2), random = ~ obstime | ID, data = aids)
# Cox proportional hazards model fit with a time-dependent covariate
fitCOX2 <- coxph(Surv(start, stop, event) ~ drug + as.numeric(drug) : obstime, 
                 data = aids, x = TRUE)
# joint model fit in which \code{timeVarT} must be specified
fitJT.ph4 <- jmodelTM(fitLME2, fitCOX2, aids, timeVarY = 'obstime', timeVarT = 'obstime')
summary(fitJT.ph4)
}
}

\keyword{survival}
\keyword{models}
