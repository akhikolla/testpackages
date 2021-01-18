#' @title
#' Predictive distribution of the max-stable process at target positions.
#' 
#' @author
#' Quentin Sebille
#' 
#' @description 
#' Computes the predictive distribution of \eqn{Y(\cdot)} at a set of ungauged positions \eqn{(s_1^*, ..., s_k^*)}, given data at gauged positions \eqn{(s_1, ..., s_n)}, by using the output of \cite{latent.fit} or \code{hkevp.fit}.
#'   
#' Two types of prediction are available for the HKEVP, as described in \cite{Shaby and Reich (2012)}. See details.
#' 
#' 
#' 
#'
#' @inheritParams extrapol.gev
#' 
#' @param predict.type 
#' Character string specifying the type of prediction. Must be one of "\code{kriging}" (default) or "\code{climat}". See details.
#' 
#' 
#'
#' @details 
#' The spatial prediction of \eqn{Y_t(s^*)} for a target site \eqn{s^*} and a realisation \eqn{t} of the process is described in \cite{Shaby and Reich (2012)}. This method involves a three-step procedure:
#' \enumerate{
#' \item Computation of the residual dependence process \eqn{\theta(\cdot)} at the target positions.
#' \item Computation of the conditional GEV parameters \eqn{(\mu^*,\sigma^*,\xi^*)} at the target sites. See the definition of the HKEVP in \cite{Reich and Shaby (2012)}.
#' \item Generation of \eqn{Y_t(s^*)} from an independent GEV distribution with parameters \eqn{(\mu^*,\sigma^*,\xi^*)}.
#' }
#' 
#' As sketched in \cite{Shaby and Reich (2012)}, two types of prediction are possible: the kriging-type and the climatological-type. These two types differ when the residual dependence process \eqn{\theta} is computed (first step of the prediction):
#' \itemize{
#' \item The kriging-type takes the actual value of \eqn{A} in the MCMC algorithm to compute the residual dependence process. The prediction will be the distribution of the maximum recorded at the specified targets.
#' \item The climatological-type generates \eqn{A} by sampling from the positive stable distribution with characteristic exponent \eqn{\alpha}, where \eqn{\alpha} is the actual value of the MCMC step. The prediction in climatological-type will be the distribution of what could happen in the conditions of the HKEVP dependence structure.
#' }
#'
#' Posterior distribution for each realisation \eqn{t} of the process and each target position \eqn{s^*} is represented with a sample where each element corresponds to a step of the MCMC procedure.
#'
#'
#'
#' @return
#' A three-dimensional array where:
#' \itemize{
#' \item Each row corresponds to a different realisation of the process (a block).
#' \item Each column corresponds to a target position.
#' \item Each slice corresponds to a MCMC step.}
#' 
#' 
#' @export
#' 
#' @references 
#' Reich, B. J., & Shaby, B. A. (2012). A hierarchical max-stable spatial model for extreme precipitation. The annals of applied statistics, 6(4), 1430. <DOI:10.1214/12-AOAS591>
#' 
#' Shaby, B. A., & Reich, B. J. (2012). Bayesian spatial extreme value analysis to assess the changing risk of concurrent high temperatures across large portions of European cropland. Environmetrics, 23(8), 638-648. <DOI:10.1002/env.2178>
#'
#'
#'
#' @examples
#' # Simulation of HKEVP:
#' sites <- as.matrix(expand.grid(1:3,1:3))
#' targets <- as.matrix(expand.grid(1.5:2.5,1.5:2.5))
#' all.pos <- rbind(sites, targets)
#' knots <- sites
#' loc <- all.pos[,1]*10
#' scale <- 3
#' shape <- 0
#' alpha <- .4
#' tau <- 1
#' ysim <- hkevp.rand(10, all.pos, knots, loc, scale, shape, alpha, tau)
#' yobs <- ysim[,1:9]
#' 
#' # HKEVP fit (omitting first site, used as target):
#' fit <- hkevp.fit(yobs, sites, niter = 1000)
#' 
#' # Extrapolation:
#' ypred <- hkevp.predict(fit, targets, predict.type = "kriging")
#' 
#' # Plot of the density and the true value for 4 first realizations:
#' # par(mfrow = c(2, 2))
#' # plot(density(ypred[1,1,]), main = "Target 1 / Year 1")
#' # abline(v = ysim[1,10], col = 2, lwd = 2)
#' # plot(density(ypred[2,1,]), main = "Target 1 / Year 2")
#' # abline(v = ysim[2,10], col = 2, lwd = 2)
#' # plot(density(ypred[1,2,]), main = "Target 2 / Year 1")
#' # abline(v = ysim[1,11], col = 2, lwd = 2)
#' # plot(density(ypred[2,2,]), main = "Target 2 / Year 2")
#' # abline(v = ysim[2,11], col = 2, lwd = 2)
#' 
#' 
#' 
#' 
hkevp.predict <- function(fit, targets, targets.covariates, predict.type = "kriging") {
  # Catching errors
  if (missing(targets.covariates)) targets.covariates <- cbind(1, targets)
  if (ncol(targets.covariates) != ncol(fit$spatial.covar))
    stop("Spatial covariates does not match between sites and targets!")
  
  if (fit$fit.type == "latent") {
    gev.estim <- extrapol.gev(fit, targets, targets.covariates)
    nyear <- nrow(fit$data)
    gevrand <- function(i) hkevp.rand(nyear, targets, targets, gev.estim$loc[i,], gev.estim$scale[i,], gev.estim$shape[i,], 1, 1)
    res.list <- tapply(1:fit$nstep, 1:fit$nstep, gevrand)
    
    RESULT <- array(NA, dim = c(nyear, ntargets, nstep))
    for (i in 1:length(res.list))
      RESULT[,,i] <- res.list[[i]]
    return(RESULT)
  }
  else if (fit$fit.type == "hkevp") {
    # Useful variables
    ntargets <- nrow(targets)
    nyear <- nrow(fit$data)
    nstep <- fit$nstep
    nknots <- nrow(fit$knots)
    dtk <- as.matrix(dist(rbind(targets, fit$knots)))[1:ntargets, -(1:ntargets)]
    if (ntargets == 1) dtk <- matrix(dtk, nrow = 1)
    GEV.extrapol <- extrapol.gev(fit, targets, targets.covariates) # Spatial extrapolation of GEV distribution at targets
    
    
    RESULT <- array(NA, dim = c(nyear, ntargets, nstep))
    
    # Function that generates one GEV(loc, scale, shape) simulation
    gev.rand <- function(loc, scale, shape) {
      if (shape == 0)
        loc - scale*log(rexp(1))
      else
        loc + scale*(rexp(1) ^ (-shape) - 1)/shape
    }
    
    # Predictive distribution
    if (predict.type == "kriging") {
      for (i in 1:nstep) {
        # Computation of the THETA process at targets
        omega <- exp(-dtk^2/(2*fit$tau[i]^2))
        omega <- sweep(omega, MARGIN = 1, STATS = rowSums(omega), FUN = "/")
        theta.targets <- (fit$A[,,i] %*% t(omega^(1/fit$alpha[i])))^fit$alpha[i]
        
        # Computation of conditional GEV parameters
        loc.star <- matrix(GEV.extrapol$loc[i,], nyear, ntargets, byrow = TRUE) + matrix(GEV.extrapol$scale[i,]/GEV.extrapol$shape[i,], nyear, ntargets, byrow = TRUE) * (theta.targets ^ matrix(GEV.extrapol$shape[i,], nyear, ntargets, byrow = TRUE) - 1)
        scale.star <- fit$alpha[i] * matrix(GEV.extrapol$scale[i,], nyear, ntargets, byrow = TRUE) * theta.targets ^ matrix(GEV.extrapol$shape[i,], nyear, ntargets, byrow = TRUE)
        shape.star <- fit$alpha[i] * matrix(GEV.extrapol$shape[i,], nyear, ntargets, byrow = TRUE)
        
        # Generation of Y
        for (t in 1:nyear) {
          for (s in 1:ntargets) {
            RESULT[t, s, i] <- gev.rand(loc.star[t,s], scale.star[t,s], shape.star[t,s])
          }
        }
      }
    }
    else if (predict.type == "climat") {
      for (i in 1:nstep) {
        # Generating A from the actual value of alpha
        unif.gen <- matrix(runif(nyear*nknots, 0, pi), nyear, nknots)
        expo.gen <- matrix(rexp(nyear*nknots, 1), nyear, nknots)
        A.gen <- (sin((1 - fit$alpha[i])*unif.gen) / expo.gen) ^ ((1 - fit$alpha[i])/fit$alpha[i] ) *
          sin(fit$alpha[i]*unif.gen) / (sin(unif.gen)) ^ (1/fit$alpha[i])
        
        # Computing the THETA process at targets:
        omega <- exp(-dtk^2/(2*fit$tau[i]^2))
        omega <- sweep(omega, MARGIN = 1, STATS = rowSums(omega), FUN = "/")
        theta.targets <- (A.gen %*% t(omega^(1/fit$alpha[i])))^fit$alpha[i]
        
        
        # Computation of conditional GEV parameters
        loc.star <- matrix(GEV.extrapol$loc[i,], nyear, ntargets, byrow = TRUE) + matrix(GEV.extrapol$scale[i,]/GEV.extrapol$shape[i,], nyear, ntargets, byrow = TRUE) * (theta.targets ^ matrix(GEV.extrapol$shape[i,], nyear, ntargets, byrow = TRUE) - 1)
        scale.star <- fit$alpha[i] * matrix(GEV.extrapol$scale[i,], nyear, ntargets, byrow = TRUE) * theta.targets ^ matrix(GEV.extrapol$shape[i,], nyear, ntargets, byrow = TRUE)
        shape.star <- fit$alpha[i] * matrix(GEV.extrapol$shape[i,], nyear, ntargets, byrow = TRUE)
        
        # Generation of Y
        for (t in 1:nyear) {
          for (s in 1:ntargets) {
            RESULT[t, s, i] <- gev.rand(loc.star[t,s], scale.star[t,s], shape.star[t,s])
          }
        }
      }
    }
    else stop("Incorrect predict.type entry!")
    return(RESULT)
  }
  else if (fit$fit.type == "dep-only") {
    
    # Useful variables
    ntargets <- nrow(targets)
    nyear <- nrow(fit$data)
    nstep <- fit$nstep
    nknots <- nrow(fit$knots)
    dtk <- as.matrix(dist(rbind(targets, fit$knots)))[1:ntargets, -(1:ntargets)]
    if (ntargets == 1) dtk <- matrix(dtk, nrow = 1)
    
    
    RESULT <- array(NA, dim = c(nyear, ntargets, nstep))
    
    # Function that generates one GEV(loc, scale, shape) simulation
    gev.rand <- function(loc, scale, shape) {
      if (shape == 0)
        loc - scale*log(rexp(1))
      else
        loc + scale*(rexp(1) ^ (-shape) - 1)/shape
    }
    
    # Predictive distribution
    if (predict.type == "kriging") {
      for (i in 1:nstep) {
        # Computation of the THETA process at targets
        omega <- exp(-dtk^2/(2*fit$tau[i]^2))
        omega <- sweep(omega, MARGIN = 1, STATS = rowSums(omega), FUN = "/")
        theta.targets <- (fit$A[,,i] %*% t(omega^(1/fit$alpha[i])))^fit$alpha[i]
        
        # Computation of conditional GEV parameters
        loc.star <- theta.targets
        scale.star <- fit$alpha[i] * theta.targets
        shape.star <- matrix(fit$alpha[i], nyear, ntargets)
        
        # Generation of Y
        for (t in 1:nyear) {
          for (s in 1:ntargets) {
            RESULT[t, s, i] <- gev.rand(loc.star[t,s], scale.star[t,s], shape.star[t,s])
          }
        }
      }
    }
    else if (predict.type == "climat") {
      for (i in 1:nstep) {
        # Generating A from the actual value of alpha
        unif.gen <- matrix(runif(nyear*nknots, 0, pi), nyear, nknots)
        expo.gen <- matrix(rexp(nyear*nknots, 1), nyear, nknots)
        A.gen <- (sin((1 - fit$alpha[i])*unif.gen) / expo.gen) ^ ((1 - fit$alpha[i])/fit$alpha[i] ) *
          sin(fit$alpha[i]*unif.gen) / (sin(unif.gen)) ^ (1/fit$alpha[i])
        
        # Computing the THETA process at targets:
        omega <- exp(-dtk^2/(2*fit$tau[i]^2))
        omega <- sweep(omega, MARGIN = 1, STATS = rowSums(omega), FUN = "/")
        theta.targets <- (A.gen %*% t(omega^(1/fit$alpha[i])))^fit$alpha[i]
        
        # Computation of conditional GEV parameters
        loc.star <- theta.targets
        scale.star <- fit$alpha[i] * theta.targets
        shape.star <- matrix(fit$alpha[i], nyear, ntargets)
        
        # Generation of Y
        for (t in 1:nyear) {
          for (s in 1:ntargets) {
            RESULT[t, s, i] <- gev.rand(loc.star[t,s], scale.star[t,s], shape.star[t,s])
          }
        }
      }
    }
    else stop("Incorrect predict.type entry!")
    return(RESULT)
  }
  else stop("Incorrect fit.type!")
}
