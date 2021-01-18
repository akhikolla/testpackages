
\name{confint}

\alias{confint.jmodelTM}
\alias{confint.jmodelMult}

\title{Obtain Confidence Intervals for Joint Model Parameters}

\description{
\code{confint} is a generic function which computes confidence intervals for parameters in models fitted by \code{jmodelTM()} or \code{jmodelMult()}.
}

\usage{
\method{confint}{jmodelTM}(object, parm, level = 0.95, \dots)
\method{confint}{jmodelMult}(object, parm, level = 0.95, \dots)
}

\arguments{
  \item{object}{an object inheriting from class \code{jmodelTM} or \code{jmodelMult}.}
  \item{parm}{a specification of which parameters are to be given confidence intervals. As currently implemented, always give confidence intervals for all regression coefficients.}
  \item{level}{the confidence level required.}
  \item{\dots}{additional arguments required. None is used in this method.}
}

\value{
A list consists of the following components:
  \item{infoLong}{a matrix with columns giving parameter estimates as well as their lower and upper confidence limits for the regression parameters of the longitudinal process. }
  \item{infoSurv}{a matrix with columns giving parameter estimates as well as their lower and upper confidence limits for the regression parameters of the survival process.}
  \item{level}{the confidence level used in computing the confidence limits.}
}

\author{
Cong Xu \email{helenxu1112@gmail.com}
}

\examples{
\dontrun{
fitLME <- lme(proth ~ Trt * obstime, random = ~ 1 | ID, data = liver)
fitCOX <- coxph(Surv(start, stop, event) ~ Trt, data = liver, x = TRUE)
fitJT.ph <- jmodelTM(fitLME, fitCOX, liver, timeVarY = 'obstime')

# 95\% confidence intervals for the joint model parameters
confint(fitJT.ph)
}
}

\keyword{methods}
