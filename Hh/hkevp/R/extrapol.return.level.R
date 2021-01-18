#' @title
#' Spatial extrapolation of a return level.
#' 
#' @author
#' Quentin Sebille
#' 
#' 
#' 
#' @description 
#' Predictive distribution of a T-years return level at ungauged positions (targets), given the output from the MCMC procedures \code{hkevp.fit} or \code{latent.fit}.
#' 
#' 
#' 
#' @inheritParams extrapol.gev
#'
#' @param period
#' An integer indicating the wished return period T.
#'
#' 
#' 
#'
#'
#' @return
#' A matrix of predictive sample. Each column corresponds to a target position and each row to a predictive draw.
#' 
#' 
#' 
#' @details 
#' Spatial extrapolation of the return level at target positions \eqn{(s^*_1, ..., s^*_k)} is a two-step procedure:
#' \itemize{
#' \item{Estimation of the predictive distribution for GEV parameters at \eqn{(s^*_1, ..., s^*_k)}, by using \code{{extrapol.gev}}.}
#' \item{Computation of the associated return level for each state of the predictive distribution.}}
#' 
#' 
#' 
#' @seealso 
#' \code{\link{extrapol.gev}}
#' 
#' 
#'
#' @examples
#' # Simulation of HKEVP:
#' sites <- as.matrix(expand.grid(1:3,1:3))
#' knots <- sites
#' loc <- sites[,1]*10
#' scale <- 1
#' shape <- .2
#' alpha <- .4
#' tau <- 1
#' ysim <- hkevp.rand(10, sites, knots, loc, scale, shape, alpha, tau)
#' 
#' # HKEVP fit:
#' fit <- hkevp.fit(ysim, sites, niter = 1000)
#' 
#' ## Extrapolation of the 100-years return level (may need more iterations and burn-in/nthin):
#' targets <- as.matrix(expand.grid(1.5:2.5,1.5:2.5))
#' pred.sample <- extrapol.return.level(100, fit, targets)
#' pred.mean <- apply(pred.sample, 2, mean)
#' pred.sd <- apply(pred.sample, 2, sd)
#' true <- return.level(100, targets[,1]*10, scale, shape)
#' # cbind(true, pred.mean, pred.sd)
#' 
#' 
#' 
extrapol.return.level <- function(period, fit, targets, targets.covariates) {
  period <- as.integer(period)
  if (missing(targets.covariates)) targets.covariates <- cbind(1, targets)
  gev <- extrapol.gev(fit, targets, targets.covariates)
  RL.chains <- matrix(NA, fit$nstep, nrow(targets))
  for (i in 1:nrow(targets))
    RL.chains[,i] <- return.level(period, gev$loc[,i], gev$scale[,i], gev$shape[,i])
  
  return(RL.chains)
}




