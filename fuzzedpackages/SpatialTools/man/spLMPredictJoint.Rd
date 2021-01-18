\name{spLMPredictJoint} 

\alias{spLMPredictJoint}

\title{
	Returns posterior predictive sample from spLM object
} 

\description{
  The function \code{spLMPredictJoint} collects posterior predictive samples
  for a set of new locations given a \code{spLM} object from the \code{spBayes} package.
}

\usage{
  spLMPredictJoint(sp.obj, pred.coords, pred.covars, start = 1, 
    end = nrow(sp.obj$p.theta.samples), thin = 1, verbose = TRUE, n.report = 100, 
    noisy = FALSE, method = "eigen")
} 

\arguments{
  \item{sp.obj}{An \code{spLM} object returned by the \code{spLM} function in the \code{spBayes} package.}
  \item{pred.coords}{An \eqn{np \times 2}{np x 2} matrix of \eqn{np} prediction
    location coordinates in \eqn{R^2}{R^2} (e.g., easting and northing).
    The first column is assumed to be easting coordinates and the second column
    northing coordinates.  }
  \item{pred.covars}{An \eqn{n \times p}{np x k} matrix of covariates matrix associated with the new locations.}
  \item{start}{Specifies the first sample included in the composition sampling.}
  \item{end}{Specifies the last sample included in the composition.
             The default is to use all posterior samples in \code{sp.obj}. }
  \item{thin}{A sample thinning factor.  The default of 1 considers all
              samples between \code{start} and \code{end}.  For example, if \code{thin = 10} then 1 in 10 samples are considered between \code{start} and
              \code{end}. }
      \item{verbose}{If \code{TRUE}, model specification and progress of the
    sampler is printed to the screen. Otherwise, nothing is printed to
    the screen.  }
  \item{n.report}{The interval to report sampling progress. }  
  \item{noisy}{If \code{TRUE}, then the posterior sample for the response is for the signal + error noise.  The default, \code{FALSE}, assumes the user wants the error-free process.}
  \item{method}{Method used to decompose covariance matrix.  Options are "chol", "eigen", and "svd" for the Cholesky, Eigen, and singular value decomposition approaches, respectively.}
} 

\value{
The function returns a \eqn{np \times B} matrix of posterior predictive samples, where \code{B} is the number of posterior samples.  The class is \code{jointPredicitveSample}.
}

\details{This function samples from the joint posterior predictive distribution of a Bayesian spatial linear model. Specifically, it is intended to be similar to the \code{spPredict} function in the \code{spBayes} except that it samples from the joint distribution instead of the marginal distribution.  However, it will only work for \code{spLM} objects and should have the same limitations as the \code{spLM} and \code{spPredict} functions.  Note that the \code{spRecover} function is called internally to recover the posterior samples form the posterior distribution of the spatial model.  
}
%\references{ R core development. "Writing R Extentions". 2007. } 

\author{Joshua French} 

\keyword{bayesian}
\keyword{joint}
\keyword{prediction} 

\seealso{spLM, spPredict, spRecover}

\examples{
# Set parameters
n <- 100
np <- 12
n.samples <- 10
burnin.start <- .5 * n.samples + 1
sigmasq <- 1
tausq <- 0.0
phi <- 1
cov.model <- "exponential"
n.report <- 5

# Generate coordinates
coords <- matrix(runif(2 * n), ncol = 2); 
pcoords <- as.matrix(expand.grid(seq(0, 1, len = 12), seq(0, 1, len = 12)))
  
# Construct design matrices
X <- as.matrix(cbind(1, coords))
Xp <- cbind(1, pcoords)

# Specify priors
starting <- list("phi" = phi, "sigma.sq"= sigmasq, "tau.sq" = tausq)
tuning <- list("phi"=0.1, "sigma.sq"=0.1, "tau.sq"=0.1)
priors.1 <- list("beta.Norm"=list(c(1, 2, 1), diag(100, 3)),
                     "phi.Unif"=c(0.00001, 10), "sigma.sq.IG"=c(1, 1))

# Generate data
B <- rnorm(3, c(1, 2, 1), sd = 10)
phi <- runif(1, 0, 10)
sigmasq <- 1/rgamma(1, 1, 1)
V <- simple.cov.sp(D = dist1(coords), cov.model, c(sigmasq, 1/phi), error.var = tausq, 
	smoothness = nu, finescale.var = 0)
y <- X \%*\% B + rmvnorm(1, rep(0, n), V) + rnorm(n, 0, sqrt(tausq))

# Create spLM object
library(spBayes)
m1 <- spBayes::spLM(y ~ X - 1, coords = coords, starting = starting,
	tuning = tuning, priors = priors.1, cov.model = cov.model,
	n.samples = n.samples, verbose = FALSE, n.report = n.report)

# Sample from joint posterior predictive distribution
y1 <- spLMPredictJoint(m1, pred.coords = pcoords, pred.covars = Xp, 
	start = burnin.start, verbose = FALSE, method = "chol")
}