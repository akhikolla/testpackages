\name{GBCEE}
\alias{GBCEE}
\title{
Generalized BCEE algorithm
}
\description{
A generalized double robust Bayesian model averaging approach to causal effect estimation. This function accommodates both binary and continuous exposures and outcomes. More details are available in Talbot and Beaudoin (2020). 
}
\usage{
GBCEE(X, Y, U, omega, niter = 5000, family.X = "gaussian",
 family.Y = "gaussian", X1 = 1, X0 = 0, priorX = NA, priorY = NA, maxsize = NA,
 OR = 20, truncation = c(0.01, 0.99), var.comp = "asymptotic", B = 200)
}
\arguments{
  \item{X}{A vector of observed values for the exposure.
}
  \item{Y}{A vector of observed values for the outcome.
}
  \item{U}{A matrix of observed values for the \code{M} potential confounding covariates, where each column contains observed values for a potential confounding factor. A recommended implementation is to only consider pre-exposure covariates.
}
  \item{omega}{The value of the hyperparameter omega in the BCEE's outcome model prior distribution. A recommended implementation is to take \code{omega}\code{ = }\code{sqrt(n)*c}, where \code{n} is the sample size and \code{c} is a user-supplied constant value. Simulation studies suggest that values of \code{c} between 100 and 1000 yield good results. 
}
  \item{niter}{The number of iterations in the Markov chain Monte Carlo model composition (MC^3) algorithm (Madigan et al. 1995). The default is 5000, but larger values are recommended when the number of potential confounding covariates is large.
}
  \item{family.X}{Distribution to be used for the exposure model. This should be \code{"gaussian"} if the exposure is continuous or \code{"binomial"} if the exposure is binary. The default is \code{"gaussian"}.
}
\item{family.Y}{Distribution to be used for the outcome model. This should be \code{"gaussian"} if the outcome is continuous or \code{"binomial"} if the outcome is binary. The default is \code{"gaussian"}.
}
\item{X1}{The value of \eqn{X1} for contrasts comparing \eqn{E[Y^{X1}]} to \eqn{E[Y^{X0}]}.  
}
\item{X0}{The value of \eqn{X0} for contrasts comparing \eqn{E[Y^{X1}]} to \eqn{E[Y^{X0}]}.
}
\item{priorX}{A vector of length \code{M} for the prior probability of inclusion of the potential confounding covariates in the exposure model (\eqn{P(\alpha^X)}). The default is 0.5 for all covariates.    
}
\item{priorY}{A vector of length \code{M} for the prior probability of inclusion of the potential confounding covariates in the outcome model. This vector multiplies BCEE's informative prior distribution (\eqn{P(\alpha^Y)}). The default is 0.5 for all covariates.  
}
\item{maxsize}{The maximum number of covariates that can be included in a given exposure or outcome model. The default is \code{M}, which does not constrain the models' size.
}
\item{OR}{A number specifying the maximum ratio for excluding models in Occam's window for the outcome modeling step. All outcome models whose posterior probability is more than \code{OR} times smaller than the largest posterior probability are excluded from the model averaging. The posterior mass of discarded models is redistributed on the remaining models. See Madigan & Raftery, 1994. The default is 20. 
}
\item{truncation}{A vector of length 2 indicating the smallest and largest values for the estimated propensity score (\eqn{P(X = 1|U)}). Values outside those bounds are truncated to the bounds. Some truncation can help reduce the impact of near positivity violations. The default is \code{c(0.01, 0.99)}. Currently, this argument is only used when \code{family.X = "binomial"}. 
}
\item{var.comp}{The method for computing the variance of the targeted maximum likelihood estimators in the BCEE algorithm. The possible values are \code{"asymptotic"}, for the efficient influence function based estimator, and \code{"boostrap"} for the nonparametric bootstrap estimator. The default is \code{"asymptotic"}.
}
\item{B}{The number of bootstrap samples when estimating the variance using the nonparametric bootstrap. The default is 200. 
}
}

\details{When both \code{Y} and \code{X} are continuous, \code{GBCEE} estimates \eqn{\Delta = E[Y^{x+1}] - E[Y^x]}, assuming a linear effect of \code{X} on \code{Y}. When \code{Y} is continuous and \code{X} is binary, \code{GBCEE} estimates \eqn{\Delta = E[Y^{X1}] - E[Y^{X0}]}. When \code{Y} is binary, \code{GBCEE} estimates both \eqn{\Delta = E[Y^{X1}] - E[Y^{X0}]} and \eqn{\Delta = E[Y^{X1}]/E[Y^{X0}]}, regardless of if \code{X} is continuous or binary. 

The \code{GBCEE} function first computes the exposure model's posterior distribution using a Markov chain Monte Carlo model composition (MC^3) algorithm (Madigan et al. 1995). The outcome model's posterior distribution is then computed using MC^3 (Madigan et al., 1995) as described in Section 3.4 of Talbot and Beaudoin (2020).

\code{GBCEE} assumes there are no missing values in the objects \code{X}, \code{Y} and \code{U}. The \code{na.omit} function which removes cases with missing data or an imputation package might be helpful.
}

\value{
\item{beta}{The model averaged estimate of the causal effect (\eqn{\hat{\Delta}}).} 
\item{stderr}{The estimated standard error of the causal effect estimate.}
\item{models.X}{A matrix giving the posterior distribution of the exposure model. Each row corresponds to an exposure model. Within each row, the first \code{M} elements are Booleans indicating the inclusion (1) or the exclusion (0) of each potential confounding factor. The last element gives the posterior probability of the exposure model.}
\item{models.Y}{A matrix giving the posterior distribution of the outcome model after applying the Occam's window. Each row corresponds to an outcome model. Within each row, the first \code{M} elements are Booleans indicating the inclusion (1) or the exclusion (0) of each potential confounding factor. The next elements are the corresponding causal effect estimate(s) and standard error(s). The last element gives the posterior probability of the outcome model.}
}

\references{
Madigan, D., York, J., Allard, D. (1995) \emph{Bayesian graphical models for discrete data}, International Statistical Review, 63, 215-232.

Madigan, D., Raftery, A. E. (1994) \emph{Model selection and accounting for model uncertainty in graphical models using Occam's window}, Journal of the American Statistical Association, 89 (428), 1535-1546.

Talbot, D., Beaudoin, C (2020) \emph{A generalized double robust Bayesian model averaging approach to causal effect estimation with application to the Study of Osteoporotic Fractures} arXiv:2003.11588 
}
\author{
Denis Talbot
}

\seealso{
\code{\link{na.omit}}.
}
\examples{
#Example:
#In this example, both U1 and U2 are potential confounding covariates.
#Both are generated as independent N(0,1).
#X is generated as a function of both U1 and U2 with a N(0,1) error.
#Y is generated as a function of X and U1 with a N(0,1) error.
#Thus, only U1 is a confounder.
#Since both X and Y are continuous, the causal contrast estimated
#by GBCEE is E[Y^{x+1}] - E[Y^{x}] assuming a linear trend.
#The true value of the causal effect is 1. 
#Unbiased estimation is possible when adjusting for U1 or
#adjusting for both U1 and U2.


#Generating the data
set.seed(418949); 
U1 = rnorm(200); 
U2 = rnorm(200);
X = 0.5*U1 + 1*U2 + rnorm(200);
Y = 1*X + 0.5*U1 + rnorm(200);

#Using GBCEE to estimate the causal exposure effect
#Very few iterations are necessary since there are only 2 covariates
results = GBCEE(X,Y,cbind(U1,U2), omega = 500*sqrt(200), niter = 50,
                family.X = "gaussian", family.Y = "gaussian");

#Causal effect estimate
results$beta;

#Estimated standard error
results$stderr;

#Results from individual models
results$models.Y;

#Posterior probability of inclusion of each covariate in the outcome model
colSums(results$models.Y[,1:2]*results$models.Y[,ncol(results$models.Y)]);
}
\keyword{causal}
\keyword{confounding}
\keyword{model average}
