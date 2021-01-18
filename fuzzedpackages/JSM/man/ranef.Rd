
\name{ranef}

\alias{ranef.jmodelTM}
\alias{ranef.jmodelMult}

\title{Extract Random Effects for Joint Models}

\description{
\code{ranef} is a generic function which extracts random effects from objects returned by \code{jmodelTM()} or \code{jmodelMult()}.
}

\usage{
\method{ranef}{jmodelTM}(object, \dots)
\method{ranef}{jmodelMult}(object, \dots)
}

\arguments{
  \item{object}{an object inheriting from class \code{jmodelTM} or \code{jmodelMult}.}
  \item{\dots}{additional arguments required. None is used in this method.}
}

\value{
A numeric matrix with rows denoting the subjects and columns the random effects. 
}

\author{
Cong Xu \email{helenxu1112@gmail.com}
}

\examples{
\dontrun{
fitLME <- lme(proth ~ Trt * obstime, random = ~ 1 | ID, data = liver)
fitCOX <- coxph(Surv(start, stop, event) ~ Trt, data = liver, x = TRUE)
fitJT.ph <- jmodelTM(fitLME, fitCOX, liver, timeVarY = 'obstime')

# random effect for the joint model
ranef(fitJT.ph)
}
}

\keyword{methods}
