\name{beam}
\alias{beam}

\title{Bayesian inference in large Gaussian graphical models}

\usage{
beam(X, type = "conditional", return.only = c("cor", "BF", "prob"),
verbose=TRUE, D=NULL)
}

\arguments{
  \item{X}{n by p data matrix}
  \item{type}{character. Either "marginal", "conditional" or "both". See Details.}
  \item{return.only}{character. Either "cor", "BF", "prob". See details.}
  \item{verbose}{logical. Whether information on progress should be be printed.}
  \item{D}{matrix. Prior marginal correlation matrix. Must be positive definite, well-conditioned and have unit variance.}
}

\value{
An object of class \code{\link{beam-class}}
}

\description{
This function carries out covariance and inverse covariance estimation within the Gaussian conjugate model. The scale matrix parameter of the inverse-Wishart is set to the identity (default), whereas the degree of freedom parameter is estimated by maximization of the marginal likelihood. The function also computes the Bayes factor and tail probability (p-values) to test the marginal or conditional independence between all pairs of variables.
}

\details{
The arguments \code{type} and \code{return.only} have essentially been introduced for computational and memory savings.
Using argument \code{type} the user may indicate whether the marginal dependencies ("marginal"), the conditional dependencies ("conditional") or both ("both") are to be inferred. On the other hand, the argument \code{return.only} is used to indicate whether the correlations ("cor"), Bayes factors ("BF") or tail probabilities ("prob") should be returned. Default is to return all three quantities for conditional dependencies.
}

\references{
Leday, G.G.R. and Richardson, S. (2019). Fast Bayesian inference in large Gaussian graphical models. \emph{Biometrics}.
}

\author{
Gwenael G.R. Leday and Ilaria Speranza
}

\examples{

# Load data
data(TCPAprad)

# beam
fit <- beam(X = TCPAprad, type="both") 

# Print summary
summary(fit)

# Extract matrix of marginal correlations
mcor(fit)[1:5, 1:5]

# Extract matrix of partial correlations
pcor(fit)[1:5, 1:5]

# Plot log-marginal likelihood of the Gaussian conjugate model
plotML(fit)

# Plot heatmap of marginal (upper triangle) and/or
# partial (lower triangle) correlation estimates
plotCor(fit)

}
