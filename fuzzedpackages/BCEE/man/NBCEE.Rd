\name{NBCEE}
\alias{NBCEE}
\title{
Naive BCEE Implementation
}
\description{
N-BCEE implementation of the BCEE algorithm. This function supports exposures that can be modeled with generalized linear models (e.g., binary, continuous or Poisson), but only continuous outcomes. 
}
\usage{
NBCEE(X, Y, U, omega, niter = 5000, nburn = 500, nthin = 10,
 maxmodelX = NA, maxmodelY = NA, family.X = "gaussian")
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
  \item{niter}{The number of post burn-in iterations in the Markov chain Monte Carlo model composition (MC^3) algorithm (Madigan et al. 1995), prior to applying thinning. The default is 5000.
}
  \item{nburn}{
The number of burn-in iterations (prior to applying thinning). The default is 500.
}
  \item{nthin}{
The thinning of the chain. The default is 10. 
}
  \item{maxmodelX}{
The maximum number of exposure models the algorithm can explore. See \code{maxmodelY} and the note below.
}
  \item{maxmodelY}{
The maximum number of distinct outcome models that the algorithm can explore. Choosing a smaller value can shorten computing time. However, choosing a value that is too small will cause the algorithm to crash. The default is \code{NA}; the maximum number of outcome models that can be explored is then set to the minimum of \code{niter + nburn} and \code{2^M}.  
}

  \item{family.X}{
A description of the error distribution and link function to be used in the model. This can be a character string naming a family function, a family function or the result of a call to a family function. (See \code{\link{family}} for details of family functions.) The default is \code{"gaussian"}
}
}
\details{
\code{NBCEE} assumes there are no missing values in the objects \code{X}, \code{Y} and \code{U}. The \code{na.omit} function which removes cases with missing data or an imputation package might be helpful.
}
\value{
\item{betas}{A vector containing the sampled values for the exposure effect.} 
\item{models.X}{A Boolean matrix identifying the sampled exposure models. See \code{models.Y}.}
\item{models.Y}{A Boolean matrix identifying the sampled outcome models. Each row corresponds to a sampled outcome model. Within each row, the \code{m}th element equals 1 if and only if the \code{m}th potential confounding covariate is included in the sampled outcome model (and 0 otherwise).}
}

\references{
Madigan, D., York, J., Allard, D. (1995) \emph{Bayesian graphical models for discrete data}, International Statistical Review, 63, 215-232.

Talbot, D., Lefebvre, G., Atherton, J. (2015) \emph{The Bayesian causal effect estimation algorithm}, Journal of Causal Inference, 3(2), 207-236.
}
\author{
Denis Talbot, Genevieve Lefebvre, Juli Atherton.
}
\note{
Variability of the exposure effect estimator is generally underestimated by the N-BCEE implementation of BCEE. The A-BCEE, which also happens to be faster, is thus preferred. Another option is to use N-BCEE with nonparametric bootstrap (B-BCEE) to correctly estimate variability. 

The difference in computing time between A-BCEE and N-BCEE is mostly explainable by the method used to compute the posterior distribution of the exposure model. In A-BCEE, this posterior distribution is calculated as a first step using \code{bic.glm}. In N-BCEE, the posterior distribution of the exposure model is computed inside the MC^3 algorithm.}

\seealso{
\code{\link{na.omit}}, \code{\link{ABCEE}}.
}

\examples{
# In this example, U1 and U2 are potential confounding covariates
# generated as independent N(0,1).
# X is generated as a function of both U1 and U2 with a N(0,1) error.
# Y is generated as a function of X and U1 with a N(0,1) error.
# Variable U1 is the only confounder.
# The causal effect of X on Y equals 1. 
# The exposure effect estimator (beta hat) in the outcome model  
# including U1 and U2 or including U1 only is unbiased.
# The sample size is n = 200.

# Generating the data
set.seed(418949); 
U1 = rnorm(200); 
U2 = rnorm(200);
X = 0.5*U1 + 1*U2 + rnorm(200);
Y = 1*X + 0.5*U1 + rnorm(200);

# Using NBCEE to estimate the causal exposure effect
n = 200;
omega.c = 500;
results = NBCEE(X,Y,cbind(U1,U2), omega = omega.c*sqrt(n),
 niter = 1000, nthin = 5, nburn = 20);

# The posterior mean of the exposure effect:
mean(results$betas);
# The posterior standard deviation of the exposure effect:
sd(results$betas);
# The posterior probability of inclusion of each covariate in the exposure model:
colMeans(results$models.X);
# The posterior distribution of the exposure model:
table(apply(results$models.X, 1, paste0, collapse = ""));
# The posterior probability of inclusion of each covariate in the outcome model:
colMeans(results$models.Y);
# The posterior distribution of the outcome model:
table(apply(results$models.Y, 1, paste0, collapse = ""));
}
\keyword{causal}
\keyword{confounding}
\keyword{model average}