#' @title
#' Spatial extrapolation of GEV parameters with the HKEVP
#' 
#' @author
#' Quentin Sebille
#' 
#' 
#' 
#' @description 
#' Predictive distributions of the GEV parameters at a set of ungauged sites (targets), given the output from the MCMC procedures \code{hkevp.fit} or \code{latent.fit}. See details.
#' 
#' 
#' 
#' 
#'
#' @param fit
#' Output from the \code{hkevp.fit} procedure.
#' 
#' @param targets
#' A matrix of real values giving the spatial coordinates of the ungauged positions. Each row corresponds to an ungauged position.
#' 
#' @param targets.covariates
#' A matrix of real values giving the spatial covariates of the ungauged positions. Must match with the covariates used in \code{\link{hkevp.fit}} or \code{\link{latent.fit}}.
#' 
#' 
#' 
#' 
#' 
#' @details 
#' Since the GEV parameters are modelled with latent Gaussian processes, spatial extrapolation of the marginal distributions at target positions \eqn{(s^*_1, ..., s^*_k)} is performed with simple kriging. Estimation is done at each MCMC step to produce a sample of the predictive distribution.
#' 
#'
#'
#' @return
#' A named list of three elements: \code{loc}, \code{scale}, \code{shape}. Each one is a matrix with columns corresponding to targets positions.
#' 
#' 
#' @seealso 
#' \code{\link{extrapol.return.level}}
#' 
#' 
#' 
#' 
#'
#' @examples
#' # Simulation of HKEVP:
#' sites <- as.matrix(expand.grid(1:3,1:3))
#' loc <- sites[,1]*10
#' scale <- 3
#' shape <- 0
#' alpha <- .4
#' tau <- 1
#' ysim <- hkevp.rand(10, sites, sites, loc, scale, shape, alpha, tau)
#' 
#' # HKEVP fit:
#' fit <- hkevp.fit(ysim, sites, niter = 1000)
#' 
#' ## Extrapolation:
#' targets <- matrix(1.5, 1, 2)
#' gev.targets <- extrapol.gev(fit, targets)
#' 
#' ## True vs predicted:
#' predicted <- sapply(gev.targets, median)
#' sd.predict <- sapply(gev.targets, sd)
#' true <- c(targets[,1]*10, scale, shape)
#' # cbind(true, predicted, sd.predict)
#' 
#' 
extrapol.gev <- function(fit, targets, targets.covariates) {
  # Default value of targets covariates and test for compatibility with sites
  if (missing(targets.covariates)) targets.covariates <- cbind(1, targets)
  if (ncol(targets.covariates) != ncol(fit$spatial.covar)) stop("Spatial covariates does not match between sites and targets!")
  all.covar <- rbind(fit$spatial.covar, targets.covariates)
  
  
  # Useful parameters
  ntargets <- nrow(targets)
  all.coord <- rbind(fit$sites, targets)
  H <- as.matrix(dist(all.coord))
  GEV.vary <- fit$spatial$vary
  nstep <- fit$nstep
  
  
  # Initializing the results
  TABLE <- matrix(NA, nstep, ntargets)
  result <- list(loc = TABLE, scale = TABLE, shape = TABLE)
  
  # Covariance function depending on the family
  if (fit$correlation == "gauss") covariance.fun <- function(h, sill, range) sill * exp(-1/2*(h/range) ^ 2)
  if (fit$correlation == "expo") covariance.fun <- function(h, sill, range) sill * exp(-h / range)
  if (fit$correlation == "mat32") covariance.fun <- function(h, sill, range) sill * (1 + sqrt(3)*h/range) * exp(-sqrt(3)*h/range)
  if (fit$correlation == "mat52") covariance.fun <- function(h, sill, range) sill * (1 + sqrt(5)*h/range  + (5/3)*(h/range) ^ 2) * exp(-sqrt(5)*h/range)
  
  
  # Simple kriging function
  simple.kriging <- function(obs, mean, covar.mat) {
    n.obs <- length(obs)
    as.vector(mean[-(1:n.obs)] + covar.mat[-(1:n.obs),1:n.obs] %*%
                solve(covar.mat[1:n.obs,1:n.obs]) %*%
                (obs - mean[1:n.obs]))
  }
  
  
  # Kriging the GEV parameters at each MCMC state after burn-in
  for (iter in 1:nstep) {
    
    # Covariance functions for the GEV parameters
    sills <- fit$spatial$sills[iter,]
    ranges <- fit$spatial$ranges[iter,]
    loc.cov <- covariance.fun(h = H, sill = sills[1], range = ranges[1])
    scale.cov <- covariance.fun(h = H, sill = sills[2], range = ranges[2])
    shape.cov <- covariance.fun(h = H, sill = sills[3], range = ranges[3])
    
    # Mean of the GEV parameters
    loc.mean <- all.covar %*% fit$spatial$beta[iter,,1]
    scale.mean <- all.covar %*% fit$spatial$beta[iter,,2]
    shape.mean <- all.covar %*% fit$spatial$beta[iter,,3]
    
    # Kriging procedure
    loc.krig <- simple.kriging(obs = fit$GEV[,1,iter], mean = loc.mean, covar.mat = loc.cov)
    if (fit$log.scale) scale.krig <- simple.kriging(obs = log(fit$GEV[,2,iter]), mean = scale.mean, covar.mat = scale.cov)
    else scale.krig <- simple.kriging(obs = fit$GEV[,2,iter], mean = scale.mean, covar.mat = scale.cov)
    shape.krig <- simple.kriging(obs = fit$GEV[,3,iter], mean = shape.mean, covar.mat = shape.cov)
    if (!GEV.vary[1]) loc.krig <- rep(fit$GEV[1,1,iter], ntargets)
    if (!GEV.vary[2]) scale.krig <- rep(log(fit$GEV[1,2,iter]), ntargets)
    if (!GEV.vary[3]) shape.krig <- rep(fit$GEV[1,3,iter], ntargets)
    
    # Saving into result table
    result$loc[iter,] <- loc.krig
    result$scale[iter,] <- scale.krig
    result$shape[iter,] <- shape.krig
    
  }
  
  if (fit$log.scale) result$scale <- exp(result$scale)
  names(result)[2] <- "scale"
  
  
  
  return(result)
}