
\name{fitted}

\alias{fitted.jmodelTM}
\alias{fitted.jmodelMult}

\title{Extract Fitted Values for Joint Models}

\description{
\code{fitted} is a generic function which extracts fitted values from objects returned by \code{jmodelTM()} or \code{jmodelMult()}.
}

\usage{
\method{fitted}{jmodelTM}(object, process = c("Longitudinal", "Survival"), 
        type = c("Marginal", "Conditional"), \dots)
\method{fitted}{jmodelMult}(object, process = c("Longitudinal", "Survival"), 
        type = c("Marginal", "Conditional"), \dots)
}

\arguments{
  \item{object}{an object inheriting from class \code{jmodelTM} or \code{jmodelMult}.}
  \item{process}{for which process the fitted values are calculated, i.e. the longitudinal or the survival process.}
  \item{type}{what type of fitted values to calculate for each process. See \bold{Details}.}
  \item{\dots}{additional arguments required. None is used in this method.}
}

\details{
We have implemented the fitted value calculation for \code{process = "Longitudinal"} but not for \code{process = "Survival"} yet as they are not well defined under the joint modeling setting. There are two types of fitted values depending on whether to compute the values conditional on the random effects. With \code{type = "Marginal"}, the fitted values are \eqn{\mathbf{X}_i^\top(t)\boldsymbol\beta} for objects returned by \code{jmodelTM()} and \eqn{\mathbf{B}^\top(t)\boldsymbol\gamma} for objects returned by \code{jmodelMult()}. With \code{type = "Conditional"}, the fitted values are \eqn{\mathbf{X}_i^\top(t)\boldsymbol\beta + \mathbf{Z}_i^\top(t)\mathbf{b}_i} for objects returned by \code{jmodelTM()} and \eqn{b_i\times\mathbf{B}^\top(t)\boldsymbol\gamma} for objects returned by \code{jmodelMult()}.
}

\value{
A numeric vector of fitted values.
}

\author{
Cong Xu \email{helenxu1112@gmail.com}
}

\examples{
\dontrun{
fitLME <- lme(proth ~ Trt * obstime, random = ~ 1 | ID, data = liver)
fitCOX <- coxph(Surv(start, stop, event) ~ Trt, data = liver, x = TRUE)
fitJT.ph <- jmodelTM(fitLME, fitCOX, liver, timeVarY = 'obstime')

# fitted values for the longitudinal process
fitted(fitJT.ph, type = "Conditional")
}
}

\keyword{methods}

