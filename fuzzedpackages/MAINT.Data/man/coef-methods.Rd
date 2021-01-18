\name{coef--methods}
\docType{methods}
\alias{coef}
\alias{coef-methods}
\alias{coef,IdtNDE-method}
\alias{coef,IdtSNDE-method}
\alias{coef,IdtNandSNDE-method}

\title{Methods for function coef in Package `MAINT.Data'}
\description{S4 methods for function coef. As in the generic coef S3 `stats' method, these methods extract parameter estimates for the models fitted to Interval Data.}

\keyword{methods}
\keyword{coef}
\keyword{Interval Data}

\usage{

\S4method{coef}{IdtNDE}(object, selmodel=BestModel(object), \dots)
\S4method{coef}{IdtSNDE}(object, selmodel=BestModel(object), ParType=c("Centr", "Direct", "All"), \dots)
\S4method{coef}{IdtNandSNDE}(object, selmodel=BestModel(object),  ParType=c("Centr", "Direct", "All"), \dots)

}

\arguments{
  \item{object}{ An object representing a model fitted to interval data.}
  \item{selmodel}{Selected model from a list of candidate models saved in object.}
  \item{ParType}{Parameterization of the Skew-Normal distribution. Only used when object has class \code{\linkS4class{IdtSNDE}} or \code{\linkS4class{IdtNandSNDE}} and in this latter case when argument \dQuote{selmodel} chooses a Skew-Normal model. \cr Alternatives are \dQuote{Centr} for centred parameters, \dQuote{Direct} for direct parameters and \dQuote{All}, for both types of parameters. See Arellano-Valle and Azzalini (2008) for details.}
  \item{\dots}{Additional arguments for method functions.}
}

\value{ A list of parameter estimates. The list components depend on the model and parametriztion assumed by the model. For Gaussian models these are respectivelly mu (vector of mean estimates) and Sigma (matrix of covariance estimates). For Skew-Normal models the components are mu, Sigma and gamma1 (one vector of skewness coefficient estimates) for the centred parametrization and the vectors ksi, and alpha, and the matrix Omega for the direct parametrization. }

\references{
Arellano-Valle, R. B. and Azzalini, A. (2008): "The centred parametrization for the multivariate skew-normal distribution". Journal of Multivariate Analysis, Volume 99, Issue 7, 1362-1382.
}

\seealso{
  \code{\link{stdEr}}, \code{\link{vcov}} 
}

\examples{

# Create an Interval-Data object containing the intervals for 899 observations 
# on the temperatures by quarter in 60 Chinese meteorological stations.

ChinaT <- IData(ChinaTemp[1:8],VarNames=c("T1","T2","T3","T4"))

ChinaT_NE <- mle(ChinaT)

# Display model estimates

print(coef(ChinaT_NE))

\dontrun{

# Estimate Skew-Normal distribution  parameters by maximum likelihood  

ChinaT_SNE <- mle(ChinaT,Model="SKNormal")

# Display model estimates

print(coef(ChinaT_SNE,ParType="Centr"))
print(coef(ChinaT_SNE,ParType="Direct"))

}

 
}

