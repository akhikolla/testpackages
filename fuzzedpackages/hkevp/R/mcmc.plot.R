#' @title
#' Markov chains plotting
#' 
#' @author
#' Quentin Sebille
#' 
#' 
#' 
#' @description 
#' Plots of the resulting Markov chains obtained by the MCMC procedures \code{hkevp.fit} or \code{latent.fit}. May be used to assess graphically convergence of the chains.
#'
#'
#'
#'
#' 
#' @param fit
#' Output from the \code{hkevp.fit} procedure.
#' 
#' @param plot.spatial
#' Logical indicating if the Markov chains of the sills and ranges hyperparameters should be plotted. FALSE by default.
#' 
#' @param mfrow
#' Optional vector of two numerical values indicating the parameter of the window plotting called by the \code{plot(...)} function.
#'
#' 
#' 
#' 
#'
#' @examples
#' # Simulation of HKEVP:
#' sites <- as.matrix(expand.grid(1:3,1:3))
#' knots <- sites
#' loc <- sites[,1]*10
#' scale <- 3
#' shape <- .2
#' alpha <- .4
#' tau <- 1
#' ysim <- hkevp.rand(10, sites, knots, loc, scale, shape, alpha, tau)
#' 
#' # HKEVP fit:
#' fit <- hkevp.fit(ysim, sites, niter = 1000)
#' 
#' # Markov chains plot:
#' mcmc.plot(fit, TRUE)
#' 
#' 
#' 
#' 
mcmc.plot <- function(fit, plot.spatial, mfrow) {
  # Default value
  if (missing(plot.spatial)) plot.spatial <- FALSE
  if (missing(mfrow)) {
    mfrow.init <- par("mfrow")
    if (plot.spatial) {
      if (fit$fit.type == "hkevp") mfrow <- c(3,3)
      else if (fit$fit.type == "dep-only") mfrow <- c(1,3)
      else if (fit$fit.type == "latent") mfrow <- c(2,3)
    } else {
      if (fit$fit.type == "dep-only") mfrow <- c(1,3)
      else if (fit$fit.type == "hkevp") mfrow <- c(2,3)
      else if (fit$fit.type == "latent") mfrow <- c(2,2)
    }
  }
  
  # Boolean that differenciates the fit types
  if (fit$fit.type == "hkevp") bool <- c(TRUE, TRUE)
  else if (fit$fit.type == "dep-only") bool <- c(FALSE, TRUE)
  else if (fit$fit.type == "latent") bool <- c(TRUE, FALSE)
  else stop("Incorrect fit.type!")
  names(bool) <- c("margins", "dep")
  par(mfrow = mfrow)
  
  if (bool[1]) {
    # Colors
    nsites <- nrow(fit$sites)[1]
    COLORS <- rgb(0:(nsites - 1)/nsites,0:(nsites - 1)/nsites,0:(nsites - 1)/nsites)
    COLORS <- rainbow(nsites)
    
    # GEV loc
    if (fit$spatial$vary[1])
      matplot(t(fit$GEV[,1,]), type = 'l', lty = 1, col = COLORS, xlab = 'Iterations', ylab = '', axes = FALSE)
    else 
      plot(fit$GEV[1,1,], type = 'l', lty = 1, col = COLORS, xlab = 'Iterations', ylab = '', axes = FALSE)
    box()
    axis(1)
    axis(2, las = 2)
    title(expression("Location "*mu))
    
    # GEV scale
    if (fit$spatial$vary[2])
      matplot(t(fit$GEV[,2,]), type = 'l', lty = 1, col = COLORS, xlab = 'Iterations', ylab = '', axes = FALSE)
    else 
      plot(fit$GEV[1,2,], type = 'l', lty = 1, col = COLORS, xlab = 'Iterations', ylab = '', axes = FALSE)
    box()
    axis(1)
    axis(2, las = 2)
    title(expression("Scale "*sigma))
    
    # GEV shape
    if (fit$spatial$vary[3])
      matplot(t(fit$GEV[,3,]), type = 'l', lty = 1, col = COLORS, xlab = 'Iterations', ylab = '', axes = FALSE)
    else 
      plot(fit$GEV[1,3,], type = 'l', lty = 1, col = COLORS, xlab = 'Iterations', ylab = '', axes = FALSE)
    box()
    axis(1)
    axis(2, las = 2)
    title(expression("Shape "*xi))
  }
  if (bool[2]) {
    # Alpha
    plot(fit$alpha, type = 'l', xlab = 'Iterations', ylab = '', axes = FALSE)
    box()
    axis(1)
    axis(2, las = 2)
    title(expression(alpha))
    
    # Tau
    plot(fit$tau, type = 'l', xlab = 'Iterations', ylab = '', axes = FALSE)
    box()
    axis(1)
    axis(2, las = 2)
    title(expression(tau))
  }
  
  # Log-likelihood
  plot(fit$llik, type = 'l', xlab = 'Iterations', ylab = '', axes = FALSE)
  box()
  axis(1)
  axis(2, las = 2)
  title(expression("log-likelihood"))
  
  
  # 4/ GEV Spatial parameters (if plot.spatial is TRUE, optional)
  if (plot.spatial & bool[1]) {
    COLORS.spat <- rainbow(3)
    
    # Sills
    matplot(fit$spatial$sills, type = 'l', lty = 1, col = COLORS.spat, xlab = 'Iterations', ylab = '', axes = FALSE)
    box(); axis(1); axis(2, las = 2)
    title(expression("Sills "*delta))
    
    
    # Ranges
    matplot(fit$spatial$ranges, type = 'l', lty = 1, col = COLORS.spat, xlab = 'Iterations', ylab = '', axes = FALSE)
    box(); axis(1); axis(2, las = 2)
    title(expression("Ranges "*lambda))
  }
  
  par(mfrow = mfrow.init)
}
