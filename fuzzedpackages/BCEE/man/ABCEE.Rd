\name{ABCEE}
\alias{ABCEE}
\title{
Approximate BCEE Implementation
}
\description{
A-BCEE implementation of the BCEE algorithm. This function supports exposures that can be modeled with generalized linear models (e.g., binary, continuous or Poisson), but only continuous outcomes. 
}
\usage{
ABCEE(X, Y, U, omega, forX = NA, niter = 5000, nburn = 500, nthin = 10,
 maxmodelY = NA, OR = 20, family.X = "gaussian")
}
\arguments{
  \item{X}{A vector of observed values for the exposure.
}
  \item{Y}{A vector of observed values for the continuous outcome.
}
  \item{U}{A matrix of observed values for the \code{M} potential confounding covariates, where each column contains observed values for a potential confounding factor. A recommended implementation is to only consider pre-exposure covariates.
}
  \item{omega}{The value of the hyperparameter omega in the BCEE's outcome model prior distribution. A recommended implementation is to take \code{omega}\code{ = }\code{sqrt(n)*c}, where \code{n} is the sample size and \code{c} is a user-supplied constant value. Simulation studies suggest that values of \code{c} between 100 and 1000 yield good results. 
}
  \item{forX}{
A Boolean vector of size \code{M}, where the \code{m}th element indicates whether or not the \code{m}th potential confounding covariate should be considered in the exposure modeling step of the BCEE algorithm. The default for \code{forX} is \code{NA}, which indicates that all potential confounding covariates should be considered in the exposure modeling step. 
}
  \item{niter}{The number of post burn-in iterations in the Markov chain Monte Carlo model composition (MC^3) algorithm (Madigan et al. 1995), prior to applying thinning. The default is 5000, but users should ensure that the value is large enough so that the number of retained samples can provide good inferences. 
}
  \item{nburn}{
The number of burn-in iterations (prior to applying thinning). The default is 500, but users should ensure that the value is large enough so that convergence of the chain is reached. An example of diagnostics of convergence of the chain is provided below.
}
  \item{nthin}{
The thinning of the chain. The default is 10, but users should ensure that the value is large enough so that there is no auto-correlation between sampled values. An example of diagnostics of absence of auto-correlation is provided below.
}
  \item{maxmodelY}{
The maximum number of distinct outcome models that the algorithm can explore. Choosing a smaller value can shorten computing time. However, choosing a value that is too small will cause the algorithm to crash and return an error message. The default is \code{NA}; the maximum number of outcome models that can be explored is then set to the minimum of \code{niter + nburn} and \code{2^M}.  
}
  \item{OR}{
A number specifying the maximum ratio for excluding models in Occam's window for the exposure modeling step (see the \code{bic.glm} help file, and Madigan & Raftery, 1994). The default is 20.  
}

  \item{family.X}{
A description of the error distribution and link function to be used in the model. This can be a character string naming a family function, a family function or the result of a call to a family function. (See \code{\link{family}} for details of family functions.) The default is \code{"gaussian"}
}
}

\details{
The \code{ABCEE} function first computes the exposure model's posterior distribution using the \code{bic.glm} function if the number of covariates is smaller than 50. Otherwise, the exact procedures depends on the value of \code{family.X}. The outcome model's posterior distribution is then computed using MC^3 (Madigan et al., 1995) as described in Talbot et al. (2015).

\code{ABCEE} assumes there are no missing values in the objects \code{X}, \code{Y} and \code{U}. The \code{na.omit} function which removes cases with missing data or an imputation package might be helpful.
}

\value{
\item{betas}{A vector containing the sampled values for the exposure effect.} 
\item{models.X}{A matrix giving the posterior distribution of the exposure model. Each row corresponds to an exposure model. Within each row, the first \code{M} elements are Booleans indicating the inclusion (1) or the exclusion (0) of each potential confounding factor. The last element gives the posterior probability of the exposure model.}
\item{models.Y}{A Boolean matrix identifying the sampled outcome models. Each row corresponds to a sampled outcome model. Within each row, the \code{m}th element equals 1 if and only if the \code{m}th potential confounding covariate is included in the sampled outcome model (and 0 otherwise).}
}

\references{
Madigan, D., York, J., Allard, D. (1995) \emph{Bayesian graphical models for discrete data}, International Statistical Review, 63, 215-232.

Madigan, D., Raftery, A. E. (1994) \emph{Model selection and accounting for model uncertainty in graphical models using Occam's window}, Journal of the American Statistical Association, 89 (428), 1535-1546.

Talbot, D., Lefebvre, G., Atherton, J. (2015) \emph{The Bayesian causal effect estimation algorithm}, Journal of Causal Inference, 3(2), 207-236.
}
\author{
Denis Talbot, Yohann Chiu, Genevieve Lefebvre, Juli Atherton.
}

\seealso{
\code{\link[BMA]{bic.glm}}, \code{\link{na.omit}}, \code{\link{NBCEE}}.
}
\examples{

#Example:
#In this example, both U1 and U2 are potential confounding covariates.
#Both are generated as independent N(0,1).
#X is generated as a function of both U1 and U2 with a N(0,1) error.
#Y is generated as a function of X and U1 with a N(0,1) error.
#Thus, only U1 is a confounder.
#The causal effect of X on Y equals 1. 
#The parameter beta associated to exposure in the outcome model
#that includes U1 and the one from the full outcome model is an
#unbiased estimator of the effect of X on Y.

#Generating the data
set.seed(418949); 
U1 = rnorm(200); 
U2 = rnorm(200);
X = 0.5*U1 + 1*U2 + rnorm(200);
Y = 1*X + 0.5*U1 + rnorm(200);

#Using ABCEE to estimate the causal exposure effect
results = ABCEE(X,Y,cbind(U1,U2), omega = 500*sqrt(200), niter = 10000, nthin = 5, nburn = 500);

##Diagnostics of convergence of the chain:
plot.default(results$betas, type = "l");
lines(smooth.spline(1:length(results$beta), results$beta), col = "blue", lwd = 2);
#The plot shows no apparent trend.
#The smoothing curve confirms that there is little or no trend,
#suggesting the chain has indeed converged before burn-in iterations ended.
#Otherwise, the value of nburn should be increased.

##Diagnostics of absence of auto-correlation
acf(results$betas, main = "ACF plot");
#Most lines are within the confidence intervals' limits, which suggests
#that there is no residual auto-correlation. If there were, the value
#of nthin should be increased.

##The number of sampled values is niter/nthin = 2000, which should be
##large enough to provide good inferences for 95\% confidence intervals.

#The posterior mean of the exposure effect:
mean(results$betas);
#The posterior standard deviation of the exposure effect:
sd(results$betas);
#The posterior inclusion probability of each covariate:
colMeans(results$models.Y);
#The posterior distribution of the outcome model:
table(apply(results$models.Y, 1, paste0, collapse = ""));

}
\keyword{causal}
\keyword{confounding}
\keyword{model average}
