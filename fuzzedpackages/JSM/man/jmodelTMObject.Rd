
\name{jmodelTMObject}

\alias{jmodelTMObject}

\title{Fitted jmodelTM Object}

\description{
An object returned by the \code{jmodelTM} function, inheriting from class \code{jmodelTM} and representing a fitted joint model for survival and longitudinal data. Objects of this class have methods for the generic functions \code{AIC}, \code{BIC}, \code{logLik}, \code{print}, \code{summary}, and \code{vcov}.
}

\value{
  The following components must be included in a legitimate \code{jmodelTM} object.
  \item{coefficients}{a list with the estimated parameters. The list is consist of the following components:}
    \describe{
      \item{beta}{the vector of estimated coefficients for the fixed effects in the linear mixed-effects model.}
      \item{phi}{the vector of estimated coefficients for the covariates other than the covariate associated with the longitudinal process in the survival model.}
      \item{alpha}{the estimated coefficient for the covariate associated with the longitudinal process in the survival model.}
      \item{Ysigma}{the estimated measurement error standard deviation for the linear mixed-effects model.}
      \item{BSigma}{the estimated variance-covariance matrix of the random effects.}
      \item{lamb}{a numeric matrix with two columns: the first column contains the unique observed survival times in ascending order; the second column contains the corresponding estimated baseline hazard values.}
      }
  \item{Vcov}{the variance-covariance matrix evaluated at the estimated parameter values.}
  \item{logLik}{the log-likelihood (the joint likelihood) value.}
  \item{est.bi}{the estimated values for the random effects}
  \item{call}{a list containing an image of the \code{jmodelTM} call that produced the object.}
  \item{numIter}{the number of iterations used in the EM algorithm.}
  \item{convergence}{the convergence indicator: if \code{"failure"}, usually more iterations are required.}
  \item{control}{the value of the \code{control} argument passed to \code{jmodelTM}.}
  \item{time.SE}{the CPU time used to compute the standard error estimates, i.e. the time use to compute the variance-covariance matrix for the parameter estimates.}
  \item{N}{the total number of repeated measurements for the longitudinal outcome.}
  \item{n}{the number of sample units.}
  \item{d}{the censoring indicator: 0 denotes censored survival time; 1 denotes observed survival time.}
  \item{rho}{the transformation parameter used for the survival model.}
}

\author{
Cong Xu \email{helenxu1112@gmail.com}
Pantelis Z. Hadjipantelis \email{pantelis@ucdavis.edu}
}

\seealso{
    \code{\link{jmodelTM}}
}

\keyword{survival}
\keyword{models}
