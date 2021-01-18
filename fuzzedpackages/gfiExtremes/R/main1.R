#' Fiducial inference for the generalized Pareto model with known threshold
#' @description Runs the MCMC sampler of the fiducial distribution for the 
#'   generalized Pareto model with known threshold.
#'
#' @param X numeric vector of data
#' @param beta vector of probabilities corresponding to the quantiles to be 
#'   estimated
#' @param threshold value of the known threshold, must be smaller than the
#'   maximum of \code{X}
#' @param gamma.init starting value for \code{gamma} in the MCMC
#' @param sigma.init starting value for \code{sigma} in the MCMC
#' @param sd.gamma standard deviation for the proposed \code{gamma} in the MCMC
#' @param sd.sigma standard deviation for the proposed \code{sigma} in the MCMC
#' @param Jnumb number of subsamples that are taken from the Jacobian
#' @param iter number of iterations per chain (burnin excluded)
#' @param burnin number of the first MCMC iterations discarded
#' @param thin thinning number for the MCMC chain. (e.g. if it is 1 no iteration 
#'   is skipped)
#' @param nchains number of MCMC chains to run
#' @param nthreads number of threads to run the chains in parallel
#' @param seeds the seeds used for the MCMC sampler; one seed per chain, or 
#'   \code{NULL} to use random seeds
#'
#' @return An object of class \code{\link[coda:mcmc]{mcmc}} if \code{nchains=1}, 
#'   otherwise an object of class \code{\link[coda:mcmc.list]{mcmc.list}}.
#'   
#' @references Damian V. Wandler & Jan Hannig. 
#'   \emph{Generalized fiducial confidence intervals for extremes}.
#'   Extremes (2012) 15:67–87.
#'   <doi:10.1007/s10687-011-0127-9>
#' 
#' @export
#' @importFrom stats na.omit
#' @importFrom doParallel registerDoParallel
#' @importFrom parallel makeCluster stopCluster
#' @importFrom foreach foreach `%dopar%`
#'
#' @examples set.seed(666L)
#' X <- rgpareto(200L, mu = 10, gamma = 0.5, sigma = 1)
#' gf <- gfigpd1(
#'   X, beta = c(0.98, 0.99), threshold = 10, 
#'   iter = 2000L, nchains = 2L, nthreads = 2L
#' ) # note: 2*2000 iterations is not enough, I'm using these settings because 
#'   # of CRAN constraints (elapsed time must be < 5s)
#' summary(gf)
#' qgpareto(c(0.98, 0.99), mu = 10, gamma = 0.5, sigma = 1)
#' rejectionRate(gf)
#' HPDinterval(gf)
#' HPDinterval(joinMCMCchains(gf))
gfigpd1 <- function(
  X, beta, threshold, 
  gamma.init = NA, sigma.init = NA, sd.gamma = NA, sd.sigma = NA, 
  Jnumb = 50L, 
  iter = 10000L, burnin = 2000L, thin = 6L,
  nchains = 4L, nthreads = parallel::detectCores(), seeds = NULL) {
  
  stopifnot(thin >= 1L, nchains >= 1L, nthreads >= 1L)
  nthreads <- min(nthreads, nchains)
  
  X <- sort(na.omit(X)) # -->> so there's no need to sort in C++
  
  if(threshold >= X[length(X)]){
    stop(
      "The value of `threshold` is larger than the maximum of `X`."
    )
  }

  n <- sum(X >= threshold)
  if(n < 3L){
    stop(
      "The sample size is too small, or the threshold is too high."
    )
  }
  
  # Initialize the default values for the tuning parameters of the MCMC chain
  if(is.na(gamma.init) || is.na(sigma.init)) {
    fit <- gpdFit(X, threshold)
    if(is.na(gamma.init)) gamma.init <- fit[1L]
    if(is.na(sigma.init)) sigma.init <- fit[2L]
  }
  if(is.na(sd.gamma)) sd.gamma <- 0.3 / log(n, 20)
  if(is.na(sd.sigma)) sd.sigma <- sigma.init * 0.3 / log(n, 20)
  
  skip.number <- thin - 1L
  number.iterations <- (skip.number + 1L) * iter + burnin
  
  if(is.null(seeds)){
    seed1 <- sample.int(2000000L, 1)
    seeds <- seed1 + 2000000L * (0L:(nchains-1L))
  }else{
    if(length(seeds) != nchains){
      stop(
        "Please specify one seed per chain."
      )
    }
    seeds <- abs(as.integer(seeds))
  }
  
  params <- c("gamma", "sigma", paste0("beta", seq_along(beta)))
  
  # run the MCMC chain
  if(nchains == 1L){
    chain <- thinChain(MCMCchainArma(
      X, beta, gamma.init, sigma.init, 
      threshold, prob = mean(X >= threshold), 
      sd.gamma, sd.sigma, # to change
      number.iterations, Jnumb, seeds[1L]
    )[-(1L:burnin), ], skip.number)
    colnames(chain) <- params
  }else{
    if(nthreads == 1L){
      chains <- vector("list", nchains)
      for(k in 1L:nchains){
        chains[[k]] <- MCMCchainArma(
          X, beta, gamma.init, sigma.init, 
          threshold, prob = mean(X > threshold), 
          sd.gamma, sd.sigma,
          number.iterations, Jnumb, seeds[k]
        )[-(1L:burnin), ]
      }
    }else{
      cl <- makeCluster(nthreads)
      registerDoParallel(cl)
      chains <- foreach(
        k = 1L:nchains, .combine = list, .multicombine = TRUE, 
        .export = "MCMCchainArma"
      ) %dopar% {
        MCMCchainArma(
          X, beta, gamma.init, sigma.init, 
          threshold, prob = mean(X > threshold), 
          sd.gamma, sd.sigma,
          number.iterations, Jnumb, seeds[k]
        )[-(1L:burnin), ]
      }
      stopCluster(cl)
    }
    chains <- lapply(chains, thinChain, skip = skip.number)
    chains <- lapply(chains, `colnames<-`, value = params)
  }
  
  if(nchains == 1L){
    out <- coda::mcmc(chain, start = burnin+1L, thin = thin)
  }else{
    out <- 
      coda::mcmc.list(lapply(chains, coda::mcmc, start = burnin+1L, thin = thin))
  }
  
  attr(out, "beta") <- beta

  return(out)
}
