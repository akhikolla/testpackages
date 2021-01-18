# @name comire.gibbs.binary
#
# @title Gibbs sampler for CoMiRe model with binary response 
# 
# @description Posterior inference via Gibbs sampler for CoMiRe model with binary response
# 
# @param y numeric vector for the response.
# @param x numeric vector for the covariate relative to the dose of exposure.
# @param mcmc a list giving the MCMC parameters. It must include the following integers: \code{nb} giving the number of burn-in iterations, \code{nrep} giving the total number of iterations, \code{thin} giving the thinning interval, \code{ndisplay} giving the multiple of iterations to be displayed on screen while the algorithm is running (a message will be printed every \code{ndisplay} iterations).
# @param prior a list containing the values of the hyperparameters. 
# It must include the following values: 
# \itemize{
# \item \code{eta}, numeric vector of size \code{J} for the Dirichlet prior on the beta basis weights, 
# \item \code{a.pi0} and \code{b.pi0}, the prior parameters of the prior beta distribution for \eqn{\pi_0},
# \item \code{J}, parameter controlling the number of elements of the Ispline basis.
# }
# @param seed seed for random initialization.
# @param max.x maximum value allowed for \code{x}.
# @param min.x minimum value allowed for \code{x}.
# @param verbose logical, if \code{TRUE} a message on the status of the MCMC algorithm is printed to the console. Default is \code{TRUE}.
# 
#' @importFrom stats dbinom rbeta
#' @importFrom gtools rdirichlet

.comire.gibbs.binary <- function(y, x, mcmc, prior, seed, max.x=max(x), min.x=min(x), verbose = TRUE){
  # prior: eta, a.pi0, b.pi0, J
  
  # internal working variables
  n <- length(y)
  print_now <- c(mcmc$nb + 1:mcmc$ndisplay*(mcmc$nrep)/mcmc$ndisplay)
  J <- prior$J # ?
  # create the objects to store the MCMC samples
  pi0 <- rep(NA, mcmc$nrep+mcmc$nb) #pi0 = pigreco_0
  pi1 <- 1
  w <- matrix(NA, mcmc$nrep+mcmc$nb, length(prior$eta))
  
  # initialize each quantity
  
  ## parameters of the model
  pi0[1] <- mean(y) 
  w[1,] <- prior$eta/sum(prior$eta) # ?
  
  ## quantity of interest
  x.grid <- seq(0, max.x, length=100)
  beta_x <- matrix(NA, mcmc$nrep+mcmc$nb, length(x.grid))
  beta_g <- matrix(NA, mcmc$nrep+mcmc$nb, length(x))
  
  ## basis expansion
  knots <- seq(min(x)+1, max.x, length=prior$J-3)
  basisX <- function(x) iSpline(x, df=3, knots = knots, Boundary.knots=c(0,max.x+1), intercept = FALSE)
  
  phiX <- basisX(x)
  phi.grid <- basisX(x.grid)
  
  ## beta_i is the interpolating function evaluated at x_i
  beta_i <- as.double(phiX %*% w[1, ])
  
  ## P0 and P1: densit? osservate
  P0 <- stats::dbinom(y, 1, pi0[1])
  P1 <- stats::dbinom(y, 1, 1)
  
  # start the MCMC simulation 
  set.seed(seed)
  for(ite in 2:(mcmc$nrep+mcmc$nb))
  {
    # 0. Print the iteration
    if(verbose)
      {
      if(ite==mcmc$nb) cat("Burn in done\n")
      if(ite %in% print_now) cat(ite, "iterations over",
                               mcmc$nrep+mcmc$nb, "\n")
      }
    
    # 1. Aggiorno d_i (dove sono rispetto a Po e Pinf)
    # 1. Update d_i marginalising out b_i from
    d <- rbinom(n, 1, prob=(beta_i*P1)/((1-beta_i)*P0 + beta_i*P1))
    
    # 2. Update b_i from the multinomial 
    b <- .labelling_bb_C(w=w[ite-1,], phi=phiX, P0=P0, P1=P1)
      #sapply(1:n, labelling_b, w[ite-1, ], phi=phiX, P0=P0, P1=P1)

    # 3. Update w from the Dirichlet and obtain an updated function beta_i
    eta.post <- as.double(prior$eta + table(factor(b,levels=1:length(w[ite-1,]))))
    w[ite, ] <- as.double(gtools::rdirichlet(1, eta.post))
    beta_i <- as.numeric(phiX %*% w[ite, ])
    beta_i[beta_i>1] <- 1
    beta_i[beta_i<0] <- 0
    
    # 4. Update pi0
    n_0 <- sum(d==0)
    y_i <- sum(y[d==0])
    pi0[ite] <- as.double(stats::rbeta(1, as.double(prior$a.pi0 + y_i),
                               as.double(prior$b.pi0 + n_0 - y_i)))
    
    # 5. Update the values of the densities in the observed points
    P0 <- stats::dbinom(y, 1, pi0[ite])
    
    # 6. compute some posterior quanities of interest
    beta_x[ite, ] <- phi.grid %*% w[ite, ]
    beta_g[ite, ] <- beta_i
  }
  post.mean.beta <- colMeans(beta_x[-c(1:mcmc$nb),])
  post.mean.betag <- colMeans(beta_g[-c(1:mcmc$nb),])
  post.mean.pi0 <- mean(pi0[-c(1:mcmc$nb)])
  
  ci.beta <- apply(beta_x[-c(1:mcmc$nb),], 2, quantile,
                   probs=c(0.025, 0.975))
  ci.pi0 <- quantile(pi0[-c(1:mcmc$nb)], probs=c(0.025, 0.975))
  
  # output
  output <- list(
    post.means = list(beta=post.mean.beta,
                    betag=post.mean.betag, pi0=post.mean.pi0),
    ci = list(beta=ci.beta, pi0=ci.pi0),
    mcmc = list(beta=beta_x, betag = beta_g, pi0=pi0,
                phiX=phiX))
  output
}

######### applico la funzione ###########
#J <- 10
#prior <- list(pi0=mean(y), dirpar=rep(1, J)/J, a=27, b=360, J=J)
#mcmc <- list(nrep=5000, nb=2000, thin=5, ndisplay=4)
#fit.bin <- modello.mistura(y, x, mcmc=mcmc, prior=prior, seed=5, max.x=max(x), min.x=min(x))