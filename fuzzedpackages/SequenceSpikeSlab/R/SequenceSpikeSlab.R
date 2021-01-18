# Make Roxygen add appropriate references for Rcpp when generating NAMESPACE
#' @useDynLib SequenceSpikeSlab
#' @importFrom Rcpp sourceCpp
NULL

# Add main documentation page for the whole package
#' Fast Exact Bayesian Inference for the Sparse Normal Means Model
#'
#' The SequenceSpikeSlab package provides fast algorithms for exact Bayesian inference in the
#' sparse normal sequence model. It implements the methods of Van Erven and Szabo, 2018.
#' Special care has been taken to make the methods scale to large data sets, and to minimize
#' numerical errors (which arise in all software because floating point numbers are
#' represented with finite precision). 
#' 
#' There are two main functions:
#' \code{\link{general_sequence_model}} and \code{\link{fast_spike_slab_beta}}.
#' 
#' For more details see the help vignette:
#' \code{vignette("SequenceSpikeSlab-vignette", package="SequenceSpikeSlab")}
#' 
#' @docType package
#' @name SequenceSpikeSlab
NULL

################################################# MAIN WRAPPERS #####################################

#' Compute marginal posterior estimates
#' 
#' This function computes marginal posterior probabilities (slab probabilities) that data points have
#' non-zero mean for the general hierarchical prior in the sparse normal sequence model. The posterior
#' mean is also provided.
#'
#' The run-time is O(n^2) on n data points, which means that doubling the size of the data
#' leads to an increase in computation time by approximately a factor of 4. Data sets of
#' size n=25,000 should be feasible within approximately 30 minutes.
#'
#' @param x Vector of n data points
#' @param sigma Standard deviation of the Gaussian noise in the data. May also be set to "auto",
#'  in which case sigma is estimated using the estimateSigma function from the selectiveInference
#'  package
#' @param slab Slab distribution. Must be either "Laplace" or "Cauchy".
#' @param log_prior Vector of length n+1 containing the logarithms of the prior probabilities pi_n(s)
#'    that the number of spikes is equal to s for s=0,...,n. It is allowed to use an unnormalized
#'    prior that does
#'    not sum to 1, because adding any constant to the log-prior probabilities does not change the
#'    result. Instead of a vector, log_prior may also be set to "beta-binomial" as a
#'    short-hand for log_prior = lbeta(beta_kappa+(0:n),beta_lambda+n-(0:n))
#'    - lbeta(beta_kappa,beta_lambda) + lchoose(n,0:n).
#' @param Laplace_lambda Parameter of the Laplace slab
#' @param Cauchy_gamma Parameter of the Cauchy slab
#' @param beta_kappa Parameter of the beta-distribution in the beta-binomial prior
#' @param beta_lambda Parameter of the beta-distribution in the beta-binomial prior. Default value=n+1
#' @param show_progress Boolean that indicates whether to show a progress bar 
#'
#' @return list (postprobs, postmean, sigma), where postprobs is a vector of marginal posterior slab
#'   probabilities that \eqn{x[i]} has non-zero mean for \eqn{i=1,...,n}; postmean is a vector with
#'   the posterior mean for the \eqn{x[i]}; and sigma is the value of sigma (this may be of
#'   interest when the sigma="auto" option is used)
#' @export
#' @importFrom selectiveInference estimateSigma
#'
#' @examples
#' # Experiments similar to those of Castilo, Van der Vaart, 2012
#' 
#' # Generate data
#' n <- 500           # sample size
#' n_signal <- 25     # number of non-zero theta
#' A <- 5             # signal strength
#' theta <- c(rep(A,n_signal), rep(0,n-n_signal))
#' x <- theta + rnorm(n, sd=1)
#' 
#' # Choose slab
#' slab <- "Laplace"
#' Laplace_lambda <- 0.5
#' 
#' # Prior 1
#' kappa1 <- 0.4   # hyperparameter
#' logprior1 <- c(0,-kappa1*(1:n)*log(n*3/(1:n)))
#' res1 <- general_sequence_model(x, sigma=1,
#'                                slab=slab,
#'                                log_prior=logprior1,
#'                                Laplace_lambda=Laplace_lambda)
#' print("Prior 1: Elements with marginal posterior probability >= 0.5:")
#' print(which(res1$postprobs >= 0.5))
#' 
#' # Prior 2
#' kappa2 <- 0.8   # hyperparameter
#' logprior2 <- kappa2*lchoose(2*n-0:n,n)
#' res2 <- general_sequence_model(x, sigma=1,
#'                                slab=slab,
#'                                log_prior=logprior2,
#'                                Laplace_lambda=Laplace_lambda)
#' print("Prior 2: Elements with marginal posterior probability >= 0.5:")
#' print(which(res2$postprobs >= 0.5))
#' 
#' # Prior 3
#' beta_kappa <- 1      # hyperparameter
#' beta_lambda <- n+1   # hyperparameter
#' res3 <- general_sequence_model(x, sigma=1,
#'                                slab=slab,
#'                                log_prior="beta-binomial",
#'                                Laplace_lambda=Laplace_lambda)
#' print("Prior 3: Elements with marginal posterior probability >= 0.5:")
#' print(which(res3$postprobs >= 0.5))
#' 
#' # Plot means for all priors
#' M=max(abs(x))+1
#' plot(1:n, x, pch=20, ylim=c(-M,M), col='green', xlab="", ylab="", main="Posterior Means")
#' points(1:n, theta, pch=20, col='blue')
#' points(1:n, res1$postmean, pch=20, col='black', cex=0.6)
#' points(1:n, res2$postmean, pch=20, col='magenta', cex=0.6)
#' points(1:n, res3$postmean, pch=20, col='red', cex=0.6)
#' legend("topright", legend=c("posterior mean 1", "posterior mean 2", "posterior mean 3",
#'                             "data", "truth"),
#'        col=c("black", "magenta", "red", "green", "blue"), pch=20, cex=0.7)
general_sequence_model <- function(x, sigma=1, slab = "Laplace",
                                   log_prior="beta-binomial",
                                   Laplace_lambda = 0.5, Cauchy_gamma = 1,
                                   beta_kappa=1, beta_lambda,
                                   show_progress=TRUE) {
  
  # Instantiate derived parameters to default values
  n <- length(x)
  if (missing(beta_lambda)) beta_lambda <- n+1
  
  beta_binomial_prior <- (length(log_prior) == 1 && log_prior == "beta-binomial")
  if (beta_binomial_prior == FALSE & length(log_prior) != n+1) {
    stop("Invalid prior")
  }
  
  # Estimate sigma if required
  if (sigma == "auto") {
    if (show_progress) {
      message("Estimating sigma... (NB. This may be slow for large n.)")
    }
    sigma <- estimateSigma(x=diag(n), y=x, intercept=FALSE, standardize=FALSE)$sigmahat
  }
  
  # Compute phi and psi
  if (slab == "Laplace") {
    logphipsi <- SSS_log_phi_psi_Laplace(x, sigma, Laplace_lambda)
  } else if (slab == "Cauchy") {
    logphipsi <- SSS_log_phi_psi_Cauchy(x, sigma, Cauchy_gamma)
  } else {
    stop("Unsupported choice of slab: ", slab)
  }
  
  # Compute marginal posterior probabilities
  if (beta_binomial_prior == FALSE) {
    postprobs <- SSS_hierarchical_prior(logphipsi, log_prior, show_progress)
  } else {
    postprobs <- SSS_hierarchical_prior_binomial(logphipsi, beta_kappa, beta_lambda)
  }
  
  # Compute posterior means
  if (slab == "Laplace") {
    postmean <- SSS_postmean_Laplace(x, logphipsi$psi, postprobs, sigma, Laplace_lambda)
  } else if (slab == "Cauchy") {
    postmean <- SSS_postmean_Cauchy(x, logphipsi$psi, postprobs, sigma, Cauchy_gamma)
  }
  
  return(list(postprobs=postprobs, postmean=postmean, sigma=sigma))
}



#' Compute marginal posterior estimates for beta-spike-and-slab prior
#' 
#' Computes marginal posterior probabilities (slab probabilities) that data points have
#' non-zero mean for the spike-and-slab prior with a Beta(beta_kappa,beta_lambda) prior
#' on the mixing parameter. The posterior mean is also provided.
#' 
#' The run-time is O(m*n^(3/2)) on n data points, which means that doubling the size of the data
#' leads to an increase in computation time by approximately a factor of 2*sqrt(2)=2.8. Data sets of
#' size n=100,000 should be feasible within approximately 30 minutes.
#'
#' @param x Vector of n data points
#' @param sigma Standard deviation of the Gaussian noise in the data. May also be set to "auto",
#'  in which case sigma is estimated using the estimateSigma function from the selectiveInference
#'  package
#' @param m The number of discretization points used is proportional to m*sqrt(n). The larger m, the
#' better the approximation, but the runtime also increases linearly with m. The default m=20 usually
#' gives sufficient numerical precision.
#' @param slab Slab distribution. Must be either "Laplace" or "Cauchy".
#' @param Laplace_lambda Parameter of the Laplace slab
#' @param Cauchy_gamma Parameter of the Cauchy slab
#' @param beta_kappa Parameter of the beta-distribution
#' @param beta_lambda Parameter of the beta-distribution. Default value=n+1
#' @param show_progress Boolean that indicates whether to show a progress bar 
#'
#' @return list (postprobs, postmean, sigma), where postprobs is a vector of marginal posterior slab
#'   probabilities that \eqn{x[i]} has non-zero mean for \eqn{i=1,...,n}; postmean is a vector with
#'   the posterior mean for the \eqn{x[i]}; and sigma is the value of sigma (this may be of
#'   interest when the sigma="auto" option is used)
#' @export
#' @importFrom selectiveInference estimateSigma
#'
#' @examples
#' # Illustrate that fast_spike_slab_beta is a faster way to compute the same results as
#' # general_sequence_model on the beta-binomial prior
#' 
#' # Generate data
#' n <- 500          # sample size
#' n_signal <- 25    # number of non-zero theta
#' A <- 5            # signal strength
#' theta <- c(rep(A,n_signal), rep(0,n-n_signal))
#' x <- theta + rnorm(n, sd=1)
#' 
#' # Choose slab
#' slab <- "Cauchy"
#' Cauchy_gamma <- 1
#' 
#' cat("Running fast_spike_slab_beta (fast for very large n)...\n")
#' res_fss <- fast_spike_slab_beta(x, sigma=1, slab=slab, Cauchy_gamma=Cauchy_gamma)
#' 
#' cat("Running general_sequence_model (slower for very large n)...\n")
#' res_gsm <- general_sequence_model(x, sigma=1, slab=slab,
#'                                   log_prior="beta-binomial", Cauchy_gamma=Cauchy_gamma)
#' 
#' cat("Maximum difference in marginal posterior slab probabilities:",
#'     max(abs(res_gsm$postprobs - res_fss$postprobs)))
#' cat("\nMaximum difference in posterior means:",
#'     max(abs(res_gsm$postmean - res_fss$postmean)), "\n")
#' 
#' # Plot means
#' M=max(abs(x))+1
#' plot(1:n, x, pch=20, ylim=c(-M,M), col='green', xlab="", ylab="",
#'      main="Posterior Means (Same for Both Methods)")
#' points(1:n, theta, pch=20, col='blue')
#' points(1:n, res_gsm$postmean, pch=20, col='black', cex=0.6)
#' points(1:n, res_fss$postmean, pch=20, col='magenta', cex=0.6)
#' legend("topright", legend=c("general_sequence_model", "fast_spike_slab_beta",
#'                             "data", "truth"),
#'        col=c("black", "magenta", "green", "blue"), pch=20, cex=0.7)
fast_spike_slab_beta <- function(x, sigma=1, m=20, slab = "Laplace",
                                     Laplace_lambda = 0.5, Cauchy_gamma = 1,
                                     beta_kappa=1, beta_lambda,
                                     show_progress=TRUE) {
  
  # Instantiate derived parameters to default values
  n <- length(x)
  if (missing(beta_lambda)) beta_lambda <- n+1
  
  # Estimate sigma if required
  if (sigma == "auto") {
    if (show_progress) {
      message("Estimating sigma... (NB. This may be slow for large n.)")
    }
    sigma <- estimateSigma(x=diag(n), y=x, intercept=FALSE, standardize=FALSE)$sigmahat
  }
  
  # Compute phi and psi
  if (slab == "Laplace") {
    logphipsi <- SSS_log_phi_psi_Laplace(x, sigma, Laplace_lambda)
  } else if (slab == "Cauchy") {
    logphipsi <- SSS_log_phi_psi_Cauchy(x, sigma, Cauchy_gamma)
  } else {
    stop("Unsupported choice of slab: ", slab)
  }
  
  # Discretize Lambda=beta distribution
  dLambda <- SSS_discretize_Lambda_beta(m, n, beta_kappa, beta_lambda)
  
  # Compute marginal posterior probabilities
  postprobs <- SSS_discrete_spike_slab(logphipsi, dLambda, show_progress)
  
  # Compute posterior means
  if (slab == "Laplace") {
    postmean <- SSS_postmean_Laplace(x, logphipsi$psi, postprobs, sigma, Laplace_lambda)
  } else if (slab == "Cauchy") {
    postmean <- SSS_postmean_Cauchy(x, logphipsi$psi, postprobs, sigma, Cauchy_gamma)
  }
  
  return(list(postprobs=postprobs, postmean=postmean, sigma=sigma))
}


################################ SUPPORTING FUNCTIONS FOR SLABS #####################################

# Internal.
# Given two vectors of quantities, returns log(exp(log_a)+exp(log_b)), but numerically stable
logsum <- function(log_a, log_b) {
  maxs <- pmax(log_a, log_b)
  delta <- exp(-abs(log_a-log_b))
  return(maxs + log1p(delta))
}

# Internal.
# f(x;s2) = 1/sqrt(2pi s2) e^-(x^2/2*s2)
# log(...) = -0.5 * (log(2pi*s2)+x^2/s2)
calculate_log_phi <- function(x, sigma) {
  return((-log(2*pi*sigma^2)-(x/sigma)^2)/2)
  # should equal dnorm(x, sd=sigma, log=TRUE)
}

#' Calculate log of phi and psi marginal densities for Laplace(lambda) slab
#'
#' Calculate log of densities phi and psi for data vector x, where
#' \deqn{phi[i] = Normal(x[i], sigma^2)}
#' \deqn{psi[i]) = E_Laplace(\theta)[Normal(x[i] - \theta, sigma^2)]}
#' 
#' @param x data vector
#' @param sigma standard deviation of observations
#' @param lambda parameter of Laplace slab density
#' 
#' @return list (phi, psi), containing logs of phi and psi densities
#' @export
#' @importFrom stats pnorm
SSS_log_phi_psi_Laplace <- function(x, sigma, lambda) {
  x_over_sigma <- x/sigma
  log_p1 <- pnorm(x_over_sigma-lambda*sigma, 0,1, log.p = TRUE)
  log_p2 <- pnorm(-x_over_sigma-lambda*sigma,0,1, log.p = TRUE)
  log_psi <- log(lambda/2)+logsum(lambda*(sigma^2*lambda-2*x)/2+log_p1, lambda*(sigma^2*lambda+2*x)/2+log_p2)
  
  return(list(phi=calculate_log_phi(x, sigma), psi=log_psi))
}


#' Calculate log of phi and psi marginal densities for Cauchy(gamma) slab
#'
#' Calculate log of densities phi and psi for data vector x, where
#' \deqn{phi[i] = Normal(x[i], sigma^2)}
#' \deqn{psi[i]) = E_Cauchy(\theta)[Normal(x[i] - \theta, sigma^2)]}
#' 
#' @param x data vector
#' @param sigma standard deviation of observations
#' @param gamma parameter of Cauchy slab density
#' 
#' @return list (phi, psi), containing logs of phi and psi densities
#' @export
#' @importFrom stats integrate
SSS_log_phi_psi_Cauchy <- function(x, sigma, gamma) {
  n <- length(x)
  log_psi <- numeric(n)
  
  mf <- 1/(2*sigma^2)
  offset <- -log(sigma*pi*gamma*sqrt(2*pi))
  
  # due to really weird scoping rules in R, the variable x_i can be defined later
  f <- function(t) {  exp(-(x_i-t)^2*mf) / (1+(t/gamma)^2) }
  
  for (i in 1:n){
    x_i<-x[i]
    res <- integrate(f, lower = -Inf, upper = Inf)
    log_psi[i] <- offset+log(res$value)
  }
  
  return(list(phi=calculate_log_phi(x, sigma), psi=log_psi))
}

#' Compute posterior means of data points for the Laplace(lambda) slab
#'
#' @param x Data vector of length n
#' @param logpsi Vector of length n that represents a preprocessed version of the data.
#'   It should contain the logs of the psi densities of the data points, as produced 
#'   by \code{\link{SSS_log_phi_psi_Laplace}}.
#' @param postprobs Vector of marginal posterior slab probabilities that \eqn{x[i]} has
#'   non-zero mean for \eqn{i=1,...,n}.
#' @param sigma standard deviation of observations
#' @param lambda parameter of Laplace slab density
#'
#' @return Vector of n posterior means
#' @export
#' @importFrom stats integrate
SSS_postmean_Laplace <- function(x, logpsi, postprobs, sigma, lambda) {
  # NB We implement this by numerical integration instead of calculating the integral analytically
  # in terms of the cdf of a Gaussian, because we believe the numerical integration to be numerically
  # more stable
  n <- length(x)
  postmean <- numeric(n)
  
  constant_factor <- lambda/(sqrt(2*pi)*sigma * 2)
  # due to really weird scoping rules in R, the variable x_i can be defined later
  f <- function(t) { t*exp(-(x_i-t)^2/(2*sigma^2) - lambda*abs(t)) }
  
  for (i in 1:n){
    x_i<-x[i]
    zeta_i <- integrate(f, lower = -Inf, upper = Inf)$value * constant_factor
    postmean[i] <- postprobs[i] * zeta_i / exp(logpsi[i])
  }
  
  return(postmean)
}


#' Compute posterior means of data points for the Cauchy(gamma) slab
#'
#' @param x Data vector of length n
#' @param logpsi Vector of length n that represents a preprocessed version of the data.
#'   It should contain the logs of the psi densities of the data points, as produced 
#'   by \code{\link{SSS_log_phi_psi_Cauchy}}.
#' @param postprobs Vector of marginal posterior slab probabilities that \eqn{x[i]} has
#'   non-zero mean for \eqn{i=1,...,n}.
#' @param sigma standard deviation of observations
#' @param gamma parameter of Cauchy slab density
#'
#' @return Vector of n posterior means
#' @export
#' @importFrom stats integrate
SSS_postmean_Cauchy <- function(x, logpsi, postprobs, sigma, gamma) {
  n <- length(x)
  postmean <- numeric(n)
  
  constant_factor <- 1/(sqrt(2*pi)*sigma * pi*gamma)
  # due to really weird scoping rules in R, the variable x_i can be defined later
  f <- function(t) { t*exp(-(x_i-t)^2/(2*sigma^2))*(1+(t/gamma)^2)^(-1) }
  
  
  for (i in 1:n){
    x_i<-x[i]
    zeta_i <- integrate(f, lower = -Inf, upper = Inf)$value * constant_factor
    postmean[i] <- postprobs[i] * zeta_i / exp(logpsi[i])
  }
  
  return(postmean)
}


################################ HIERARCHICAL ########################################

#' Compute marginal posterior probabilities (slab probabilities) that data points have
#'   non-zero mean for the hierarchical prior.
#'
#' @param log_phi_psi List \{logphi, logpsi\} containing two vectors of the same length n
#'   that represent a preprocessed version of the data. logphi and logpsi should contain
#'   the logs of the phi and psi densities of the data points, as produced for instance
#'   by \code{\link{SSS_log_phi_psi_Laplace}} or \code{\link{SSS_log_phi_psi_Cauchy}}
#' @param logprior vector of length n+1 with components logprior[p]=log(pi_n(p)) for
#'   \eqn{p=0,...,n}
#' @param show_progress Boolean that indicates whether to show a progress bar
#'
#' @return Returns a vector with marginal posterior slab probabilities that \eqn{x[i]} has
#'   non-zero mean for \eqn{i=1,...,n}.
#' @export
SSS_hierarchical_prior <- function(log_phi_psi, logprior, show_progress = TRUE) {
  return(HierarchicalPriorC(log_phi_psi$phi, log_phi_psi$psi, logprior, show_progress))
}

#' Compute marginal posterior probabilities (slab probabilities) that data points have
#'   non-zero mean using the general hierarchical prior algorithm, but specialized to
#'   the Beta[kappa,lambda]-binomial prior. This function is equivalent to calling
#'   \code{\link{SSS_hierarchical_prior}} with
#'   logprior = lbeta(kappa+(0:n),lambda+n-(0:n)) - lbeta(kappa,lambda) + lchoose(n,0:n),
#'   but more convenient when using the Beta[kappa,lambda]-binomial prior and with a
#'   minor interior optimization that avoids calculating the choose explicitly.
#'
#' @param log_phi_psi List \{logphi, logpsi\} containing two vectors of the same length n
#'   that represent a preprocessed version of the data. logphi and logpsi should contain
#'   the logs of the phi and psi densities of the data points, as produced for instance
#'   by \code{\link{SSS_log_phi_psi_Laplace}} or \code{\link{SSS_log_phi_psi_Cauchy}}
#' @param kappa First parameter of the beta-distribution
#' @param lambda Second parameter of the beta-distribution
#' @param show_progress Boolean that indicates whether to show a progress bar
#'
#' @return Returns a vector with marginal posterior slab probabilities that \eqn{x[i]} has
#'   non-zero mean for \eqn{i=1,...,n}.
#' @export
SSS_hierarchical_prior_binomial <- function(log_phi_psi, kappa, lambda, show_progress = TRUE) {
  n <- length(log_phi_psi$phi)
  # logprior without the lchoose part
  logprior <- lbeta(kappa+(0:n),lambda+n-(0:n)) - lbeta(kappa,lambda)
  # Then tell HierarchicalPriorC to skip subtracting the lchoose part
  return(HierarchicalPriorC(log_phi_psi$phi, log_phi_psi$psi, logprior, show_progress, divideByBinom=FALSE))
}

################################ DISCRETIZED ########################################


#' Creates a vector of uniformly spaced grid points in the beta parametrization
#' Ensures the number of generated grid points is >= mingridpoints (which does
#' not have to be integer), and that their number is always odd so there is always 
#' a grid point at pi/4.
#'
#' @param minngridpoints Minimum number of grid points 
#'
#' @return Vector of betagrid points
#' @export
SSS_make_beta_grid <- function(minngridpoints) {
  ngridpoints <- ceiling(minngridpoints)
  if (ngridpoints %% 2 == 0) ngridpoints <- ngridpoints + 1
  step <- pi/2/ngridpoints
  return(seq(from=step/2,to=pi/2,by=step)) # the values of the discrete beta
}

#' Given a prior Lambda on the alpha-parameter in the spike-and-slab model,
#' make a discretized version of Lambda that is only supported on a grid of
#' approximately m * sqrt(n) discrete values of alpha. This discretized
#' version of Lambda is required as input for
#' \code{\link{SSS_discrete_spike_slab}}.
#' NB Lambda needs to satisfy a technical condition from the
#' paper that guarantees its density does not vary too rapidly. For
#' Lambda=Beta(kappa,lambda) use \code{\link{SSS_discretize_Lambda_beta}}
#' instead.
#'
#' @param m A multiplier for the number of discretization points
#' @param n The sample size
#' @param log_Lambda_cdf A function that takes as input a value of alpha and
#' calculates the log of the cumulative distribution function of Lambda at alpha
#'
#' @return List (alpha_grid, log_probs), where alpha_grid is a vector with the
#' generated grid points, and log_probs are the logs of the prior probabilities
#' of these grid points for the discretized Lambda prior.
#' @export 
SSS_discretize_Lambda <- function(m=20, n, log_Lambda_cdf) {
  f <- function(x0, x1) {
    l0 <- log_Lambda_cdf(x0)
    l1 <- log_Lambda_cdf(x1)
    return(log1p(-exp(l0-l1))+l1)
  }
  discretize_Lambda_internal(m, n, f)
}

# Internal
# m is a multiplier for the number of discretization points
# n is the sample size.
# number of grid points will be proportional to m*sqrt(n)
# log_Lambda_cdf_delta(x0,x1) is log(Lambda(x1)-Lambda(x0)) where Lambda
# is the cdf of the Lambda prior on alpha.
discretize_Lambda_internal <- function(m, n, log_Lambda_cdf_delta) {
  beta_grid <- SSS_make_beta_grid(2*(m+1)*ceiling(sqrt(n))+1)
  n_g = length(beta_grid)
  log_alpha_prior <- numeric(n_g)
  a0 <- 0
  for (i in 1:n_g) {
    a1 <- ifelse(i < n_g, sin(0.5*(beta_grid[i+1] + beta_grid[i]))^2, 1)
    log_alpha_prior[i] <- log_Lambda_cdf_delta(a0, a1)
    a0 <- a1
  }
  alpha_grid <- sin(beta_grid)^2
  return(list(alpha_grid = alpha_grid, log_probs = log_alpha_prior))
}


#' Given prior Lambda=Beta(kappa,lambda) on the alpha-parameter in the
#' spike-and-slab model, make a discretized version of Lambda that is only
#' supported on a grid of approximately m * sqrt(n) discrete values of alpha.
#' This discretized version of Lambda is required as input for
#' SSS_discrete_spike_slab.
#'
#' @param m A multiplier for the number of discretization points
#' @param n The sample size
#' @param kappa Parameter of the prior. Needs to be at least 0.5.
#' @param lambda Parameter of the prior. Needs to be at least 0.5.
#'
#' @return List (alpha_grid, log_probs), where alpha_grid is a vector with the
#' generated grid points, and log_probs are the logs of the prior probabilities
#' of these grid points for the discretized Lambda prior.
#' @export
SSS_discretize_Lambda_beta <- function(m=20, n, kappa, lambda) {
  # [This works by constructing a grid for sample size n' = n + (kappa-1/2+lambda-1/2) according to the description
  # in the paper. Then initializing prior to beta(1/2,1/2) prior (which is uniform) and fast-forwarding it
  # to the posterior after observing kappa-1/2 fake ones and lambda-1/2 fake zeros.]
  if (kappa < 0.5) { stop("kappa needs to be at least 0.5") }
  if (lambda < 0.5) { stop("lambda needs to be at least 0.5") }
  # Normal alpha grid
  effective_n <- n + kappa + lambda - 1
  k <- 2*(m + 1)*ceiling(sqrt(effective_n)) + 1
  alpha_grid <- sin(SSS_make_beta_grid(k))^2
  k <- length(alpha_grid)
  # Start with uniform prior in beta-parametrization
  log_alpha_prior <- rep(-log(k),k)
  # Now add kappa-0.5 ones and lambda-0.5 fake zeros and compute posterior, which will be our real prior
  for (j in 1:k) {
    log_alpha_prior[j] <- log_alpha_prior[j] + (kappa-0.5)*log(alpha_grid[j]) + (lambda-0.5)*log(1-alpha_grid[j])
  }
  # Normalize
  lognorm <- -Inf
  for (j in 1:k) {
    if (lognorm < log_alpha_prior[j]) {
      lognorm <- log1p(exp(lognorm-log_alpha_prior[j]))+log_alpha_prior[j]
    } else {
      lognorm <- log1p(exp(log_alpha_prior[j]-lognorm))+lognorm
    }
  }
  log_alpha_prior <- log_alpha_prior - lognorm
  
  return(list(alpha_grid = alpha_grid, log_probs = log_alpha_prior))
}

#' Compute marginal posterior probabilities (slab probabilities) that data points have
#'   non-zero mean for the discretized spike-and-slab prior.
#'
#' @param log_phi_psi List \{logphi, logpsi\} containing two vectors of the same length n
#'   that represent a preprocessed version of the data. logphi and logpsi should contain
#'   the logs of the phi and psi densities of the data points, as produced for instance
#'   by \code{\link{SSS_log_phi_psi_Laplace}} or \code{\link{SSS_log_phi_psi_Cauchy}}
#' @param dLambda Discretized Lambda prior, as generated by either
#' discretize_Lambda or discretize_Lambda_beta.
#' @param show_progress Boolean that indicates whether to show a progress bar
#'
#' @return Returns a vector with marginal posterior slab probabilities that \eqn{x[i]} has
#'   non-zero mean for \eqn{i=1,...,n}.
#' @export
SSS_discrete_spike_slab <- function(log_phi_psi, dLambda, show_progress = TRUE) {
  return(DiscreteSpikeSlabPriorC(log_phi_psi$phi, log_phi_psi$psi, dLambda$alpha_grid, dLambda$log_probs, show_progress))
}
