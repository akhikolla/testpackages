
\name{residuals}

\alias{residuals.jmodelTM}
\alias{residuals.jmodelMult}

\title{Extract Residuals for Joint Models}

\description{
\code{residuals} is a generic function which extracts residuals from objects returned by \code{jmodelTM()} or \code{jmodelMult()}.
}

\usage{
\method{residuals}{jmodelTM}(object, process = c("Longitudinal", "Survival"), 
            type = c("Marginal", "Conditional", "Standardized-Marginal", 
            "Standardized-Conditional"), \dots)
\method{residuals}{jmodelMult}(object, process = c("Longitudinal", "Survival"), 
            type = c("Marginal", "Conditional", "Standardized-Marginal", 
            "Standardized-Conditional"), \dots)
}

\arguments{
  \item{object}{an object inheriting from class \code{jmodelTM} or \code{jmodelMult}.}
  \item{process}{for which process the residuals are calculated, i.e. the longitudinal or the survival process.}
  \item{type}{what type of residuals to calculate for each process. See \bold{Details}.}
  \item{\dots}{additional arguments required. None is used in this method.}
}

\details{
We have implemented the residual calculation for \code{process = "Longitudinal"} but not for \code{process = "Survival"} yet as they are not well defined under the joint modeling setting. There are four types of residuals depending on whether to compute the values conditional on the random effects and whether to standardize the residuals. Please refer to Nobre and Single (2007) for details. 

With \code{type = "Marginal"}, the residuals are \eqn{\varepsilon_{ij} = Y_{ij} - \mathbf{X}_{ij}^\top\boldsymbol\beta} for objects returned by \code{jmodelTM()} and \eqn{\varepsilon_{ij} = Y_{ij} - \mathbf{B}^\top(t_{ij})\boldsymbol\gamma} for objects returned by \code{jmodelMult()}. With \code{type = "Conditional"}, the residuals are \eqn{\varepsilon_{ij} = Y_{ij} - \mathbf{X}_{ij}^\top\boldsymbol\beta - \mathbf{Z}_{ij}^\top\mathbf{b}_i} for objects returned by \code{jmodelTM()} and \eqn{\varepsilon_{ij} = Y_{ij} - b_i\times\mathbf{B}^\top(t_{ij})\boldsymbol\gamma} for objects returned by \code{jmodelMult()}. If \code{type = "Standardized-Marginal"} or \code{type = "Standardized-Conditional"}, the above defined residuals are divided by the estimated standard deviation of the corresponding error term.
}

\value{
A numerc vector of residual values.
}

\references{
Nobre, J. S. and Singer, J. M. (2007) Residuals analysis for linear mixed models. \emph{Biometrical Jounral} \bold{49}(6), 863--875.
}

\author{
Cong Xu \email{helenxu1112@gmail.com}
}

\seealso{
\code{\link{fitted.jmodelTM}}, 
\code{\link{fitted.jmodelMult}}
}

\examples{
\dontrun{
fitLME <- lme(proth ~ Trt * obstime, random = ~ 1 | ID, data = liver)
fitCOX <- coxph(Surv(start, stop, event) ~ Trt, data = liver, x = TRUE)
fitJT.ph <- jmodelTM(fitLME, fitCOX, liver, timeVarY = 'obstime')

# residuals for the longitudinal process
residuals(fitJT.ph, type = "Standardized-Conditional")
}
}

\keyword{methods}
