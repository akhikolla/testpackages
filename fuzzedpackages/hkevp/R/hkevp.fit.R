#' @title
#' Fitting procedure of the HKEVP with MCMC algorithm
#' 
#' @author
#' Quentin Sebille
#' 
#' 
#' 
#' @description 
#' Metropolis-within-Gibbs algorithm that returns samples from posterior distribution of all the parameters of the HKEVP.
#' 
#' Most of the input parameters have default values, so that the procedure can be easily handled. However, convergence of the Markov chains should be assessed by using \code{\link{mcmc.plot}} for instance. The experimented user can set initial states, prior hyperparameters along with the magnitude of the MCMC jumps.
#' 
#' 
#' 
#' @inheritParams latent.fit
#' 
#' @param knots
#' The coordinates of the knots in the HKEVP. By default, the positions of the knots coincide with the positions of the sites.
#' 
#' @param fit.margins
#' A logical that indicates if the GEV parameters should be fitted along with the dependence structure. TRUE by default.
#' 
#' @param mcmc.init 
#' A named list indicating the initial states of the chains. See details.
#' 
#' @param mcmc.prior 
#' A named list indicating the hyperparameters of the prior distributions. See details.
#' 
#' @param mcmc.jumps 
#' A named list indicating the amplitude of the jumps to propose the MCMC candidates. See details.
#' 
#' 
#' 
#' 
#' 
#' @details 
#' Details of the MCMC procedure are presented in \cite{Reich and Shaby (2012)}. This function follows the indications and the choices of the authors, with the exception of several small changes:
#' \itemize{
#' \item{The scale parameter \eqn{\sigma} can be modelled like the two other marginal parameters as in \cite{Davison et al. (2012)} or by its logarithm as in \cite{Reich and Shaby (2012)}. For this, use the argument \code{log.scale}, set to FALSE by default.}
#' 
#' \item{The Inverse-Gamma prior distributions defined for the bandwith parameter \eqn{\tau} and for the ranges \eqn{\lambda} of the latent processes are replaced by a Beta distribution over the interval \eqn{[0,2D_{max}]}, where \eqn{D_{max}} stands for the maximum distance between two sites.}}
#' 
#' The procedure can be used normally with \code{fit.margins = TRUE} (default) or by assuming that the observed process had GEV(1,1,1) margins already and thus ignoring the marginal estimation.
#' 
#' If the margins are estimated and the parameters are assumed spatially-varying, the user can provide spatial covariates to fit the mean of the latent Gaussian processes. Recall for instance for the GEV location parameter that:
#' \deqn{\mu(s) = \beta_{0,\mu} + \beta_{1,\mu} c_1(s) + ... + \beta_{p,\mu} c_p(s) ~.}
#' The given matrix \code{spatial.covariates} that represents the \eqn{c_i(s)} elements should have the first column filled with ones to account for the intercept \eqn{\beta_0}.
#' 
#' The arguments \code{mcmc.init}, \code{mcmc.prior} and \code{mcmc.jumps} are named list that have default values. The user can make point changes in these arguments, by setting \code{mcmc.init = list(alpha = .5)} for instance, but must respect the constraints of each element:
#' \itemize{
#' \item{\code{mcmc.init}. All elements are of length one. The possibilities are:
#' \itemize{
#' \item{\code{loc}, \code{scale} and \code{shape} (GEV parameters).}
#' \item{\code{range} and \code{sill} of the correlation functions.}
#' \item{\code{alpha}, \code{tau}, \code{A} and \code{B}, the dependence parameters and conditional variables of the HKEVP.}}}
#' \item{mcmc.prior. The possible elements are:
#' \itemize{
#' \item{\code{constant.gev}: a \eqn{2 \times 3} matrix of normal parameters for spatially-constant \eqn{\mu}, \eqn{\sigma} and \eqn{\xi}. The first row are the means, the second are the standard deviations.}
#' \item{\code{beta.sd}: the normal sd prior of all \eqn{\beta} parameters (a single value).}
#' \item{\code{range}, \code{alpha} and \code{tau}: the two Beta parameters.}
#' \item{\code{sill}: the two Inverse-Gamma parameters.}}}
#' \item{mcmc.jumps. The possible elements are:
#' \itemize{
#' \item{\code{gev} and \code{range}: a vector of length 3 (for each GEV parameter).}
#' \item{\code{tau}, \code{alpha}, \code{A}, \code{B}: single values for each.}}}}
#' 
#' 
#' 
#' 
#'
#'
#'
#'
#' @return
#' A named list with following elements:
#' \itemize{
#' \item{\code{GEV}: the Markov chains associated to the GEV parameters. The dimensions of the array correspond respectively to the sites positions, the three GEV parameters and the states of the Markov chains.}
#' \item{\code{alpha}: the Markov chain associated to the dependence parameter \eqn{\alpha}.}
#' \item{\code{tau}: the Markov chain associated to the dependence parameter \eqn{\tau}.}
#' \item{\code{A}: the Markov chains associated to the positive stable random effect per site and per block. The dimensions correspond respectively to the indices of blocks, the knots positions and the states of the Markov chains.}
#' \item{\code{llik}: the log-likelihood of the model for each step of the algorithm.}
#' \item{\code{time}: time (in sec) spent for the fit.}
#' \item{\code{spatial}: a named list with four elements linked to the GEV spatially-varying parameters:
#' \itemize{
#' \item{\code{vary}: the argument \code{gev.vary}.}
#' \item{\code{beta}: the \eqn{\beta} parameters for each GEV parameter. The dimensions correspond respectively to the steps of the Markov chains, the \eqn{p} spatial covariates and the GEV parameters}
#' \item{\code{sills}: the Markov chains associated to the sills in the correlation functions of the latent Gaussian processes.}
#' \item{\code{ranges}: the Markov chains associated to the ranges in the correlation functions of the latent Gaussian processes.}}}
#' \item{\code{data}: the data fitted.}
#' \item{\code{sites}: the sites where the data are observed.}
#' \item{\code{knots}: the set of knots.}
#' \item{\code{spatial.covariates}: the spatial covariates.}
#' \item{\code{correlation}: the type of correlation function for the marginal latent processes.}
#' \item{\code{nstep}: the number of steps at the end of the routine after burn-in and thinning.}
#' \item{\code{log.scale}: a boolean indicating if the scale parameter has been modelled via its logarithm.}
#' \item{\code{fit.type}: either "hkevp" or "dep-only" character string to specify the type of fit.}}
#' 
#' If \code{fit.margins} is false, only the dependence-related elements are returned.
#' 
#' 
#' 
#' @seealso latent.fit
#' 
#' 
#' @export
#'
#'
#' @references 
#' Reich, B. J., & Shaby, B. A. (2012). A hierarchical max-stable spatial model for extreme precipitation. The annals of applied statistics, 6(4), 1430. <DOI:10.1214/12-AOAS591>
#' 
#' Stephenson, A. G. (2009) High-dimensional parametric modelling of multivariate extreme events. Aust. N. Z. J Stat, 51, 77-88. <DOI:10.1111/j.1467-842X.2008.00528.x>
#' 
#' Davison, A. C., Padoan, S. A., & Ribatet, M. (2012). Statistical modeling of spatial extremes. Statistical Science, 27(2), 161-186. <DOI:10.1214/11-STS376>
#' 
#' 
#' 
#' 
#' @examples
#' # Simulation of HKEVP:
#' set.seed(1)
#' sites <- as.matrix(expand.grid(1:3,1:3))
#' loc <- sites[,1]*10
#' scale <- 3
#' shape <- 0
#' alpha <- .4
#' tau <- 1
#' ysim <- hkevp.rand(10, sites, sites, loc, scale, shape, alpha, tau)
#' 
#' # HKEVP fit:
#' fit <- latent.fit(ysim, sites, niter = 1000)
#' 
#' 
#' 
hkevp.fit <- function(y, sites, knots, niter, nburn, nthin, quiet, trace, fit.margins, gev.vary, spatial.covariates, log.scale, correlation, mcmc.init, mcmc.prior, mcmc.jumps)
{
  ##########################################
  ## Default values for missing arguments ##
  ##########################################
  # General arguments
  if (missing(knots)) knots <- sites
  if (!is.matrix(y)) stop("Argument y must be a matrix!")
  if (!is.matrix(sites)) stop("Argument sites must be a matrix!")
  if (!is.matrix(knots)) stop("Argument knots must be a matrix!")
  if (missing(nburn)) nburn <- 0
  if (missing(nthin)) nthin <- 1
  if (missing(quiet)) quiet <- FALSE
  if (missing(trace)) trace <- max(floor(niter/10), 1)
  if (sum(colSums(!is.na(y)) == 0) > 0) stop("Argument y must have at least one observation per column!")
  if (ncol(y) != nrow(sites)) stop("Argument y and sites do not match!")
  if (nburn < 0 | nburn >= niter) stop("Invalid nburn parameter!")
  if (nthin <= 0 | nthin >= niter) stop("Invalid nthin parameter!")
  if (trace <= 0) trace <- 1
  trace <- floor(trace) # Protection against possible seg.fault!
  
  # Margins
  if (missing(fit.margins)) fit.margins <- TRUE
  if (missing(gev.vary)) gev.vary <- c(TRUE, TRUE, FALSE)
  if (length(gev.vary)!=3 | typeof(gev.vary)!="logical") stop("Invalid gev.vary parameter!")
  if (missing(spatial.covariates)) spatial.covariates <- cbind(1, sites)
  if (nrow(spatial.covariates) != nrow(sites)) stop("Arguments sites and spatial.covariates do not match!")
  if (missing(log.scale)) log.scale <- FALSE
  if (missing(correlation)) correlation <- "mat32"
  if (!(correlation %in% c("expo", "gauss", "mat32", "mat52"))) stop("Argument correlation must be one of 'expo', 'gauss', 'mat32' or 'mat52'!")
  
  
  # Initial states
  loc.init <- median(y, na.rm = TRUE)
  scale.init <- 5
  shape.init <- .01
  range.init <- max(dist(sites))/2
  sill.init <- 1
  alpha.init <- .25
  tau.init <- min(dist(knots))
  A.init <- exp(2)
  B.init <- .5
  if (missing(mcmc.init)) mcmc.init <- list()
  if (!is.list(mcmc.init)) stop("Argument mcmc.init must be a list!")
  if (is.null(mcmc.init$loc)) mcmc.init$loc <- loc.init
  if (is.null(mcmc.init$scale)) mcmc.init$scale <- scale.init
  if (is.null(mcmc.init$shape)) mcmc.init$shape <- shape.init
  if (is.null(mcmc.init$range)) mcmc.init$range <- range.init
  if (is.null(mcmc.init$sill)) mcmc.init$sill <- sill.init
  if (is.null(mcmc.init$alpha)) mcmc.init$alpha <- alpha.init
  if (is.null(mcmc.init$tau)) mcmc.init$tau <- tau.init
  if (is.null(mcmc.init$A)) mcmc.init$A <- A.init
  if (is.null(mcmc.init$B)) mcmc.init$B <- B.init
  if (any(lapply(mcmc.init, length) != 1)) stop("Argument mcmc.init is ill-defined!")
  if (mcmc.init$shape == 0) mcmc.init$shape <- .001  # Cannot account for exact 0 shape
  if (mcmc.init$tau<=0 | mcmc.init$tau>max(dist(sites))*2)
    stop("Invalid initial state for tau: should be between 0 and 2*Dmax!")
  if (mcmc.init$range<=0 | mcmc.init$range>max(dist(sites))*2)
    stop("Invalid initial state for range: should be between 0 and 2*Dmax!")
  if (mcmc.init$sill <= 0) stop("Invalid initial state for sills: must be positive!")
  if (mcmc.init$alpha<=0 |mcmc.init$alpha > 1) stop("Invalid initial state for alpha: must be in (0,1]!")
  if (log.scale)  mcmc.init$scale <- log(mcmc.init$scale)
  
  # Prior distributions
  constant.gev.prior <- rbind(rep(0,3), c(10, 1, .25))  # Normal
  beta.sd.prior <- 100  # Normal
  sill.prior <- c(.1,.1)  # Inverse-Gamma
  range.prior <- c(2,5)  # Beta
  alpha.prior <- c(1,1)  # Beta
  tau.prior <- c(2,5)  # Beta
  illdefined.prior <- FALSE
  if (missing(mcmc.prior)) mcmc.prior <- list()
  if (!is.list(mcmc.prior)) stop("Argument mcmc.prior must be a list!")
  if (is.null(mcmc.prior$constant.gev)) mcmc.prior$constant.gev <- constant.gev.prior
  if (is.null(mcmc.prior$beta.sd)) mcmc.prior$beta.sd <- beta.sd.prior
  if (is.null(mcmc.prior$range)) mcmc.prior$range <- range.prior
  if (is.null(mcmc.prior$sill)) mcmc.prior$sill <- sill.prior
  if (is.null(mcmc.prior$alpha)) mcmc.prior$alpha <- alpha.prior
  if (is.null(mcmc.prior$tau)) mcmc.prior$tau <- tau.prior
  if (any(dim(mcmc.prior$constant.gev) != c(2,3))) illdefined.prior <- TRUE
  if (length(mcmc.prior$beta.sd) != 1)  illdefined.prior <- TRUE
  if (length(mcmc.prior$range) != 2)  illdefined.prior <- TRUE
  if (length(mcmc.prior$sill) != 2)  illdefined.prior <- TRUE
  if (length(mcmc.prior$alpha) != 2)  illdefined.prior <- TRUE
  if (length(mcmc.prior$tau) != 2)  illdefined.prior <- TRUE
  if (illdefined.prior) stop("Argument mcmc.prior is ill-defined!")
  mcmc.prior$sill <- matrix(mcmc.prior$sill, 2, 3)
  mcmc.prior$range <- matrix(mcmc.prior$range, 2, 3)
  
    
  # Jumps length to generate candidates
  gev.jumps <- c(1, .1, .01)
  range.jumps <- c(1, 1, 1)
  tau.jumps <- .02
  alpha.jumps <- .01
  A.jumps <- 1
  B.jumps <- .25
  illdefined.jumps <- FALSE
  if (missing(mcmc.jumps)) mcmc.jumps <- list()
  if (!is.list(mcmc.jumps)) stop("Argument mcmc.jumps must be a list!")
  if (is.null(mcmc.jumps$gev)) mcmc.jumps$gev <-  gev.jumps
  if (is.null(mcmc.jumps$range)) mcmc.jumps$range <- range.jumps
  if (is.null(mcmc.jumps$tau)) mcmc.jumps$tau <- tau.jumps
  if (is.null(mcmc.jumps$alpha)) mcmc.jumps$alpha <- alpha.jumps
  if (is.null(mcmc.jumps$A)) mcmc.jumps$A <- A.jumps
  if (is.null(mcmc.jumps$B)) mcmc.jumps$B <- B.jumps
  if (length(mcmc.jumps$gev) != 3) illdefined.jumps <- TRUE
  if (length(mcmc.jumps$range) != 3) illdefined.jumps <- TRUE
  if (length(mcmc.jumps$tau) != 1) illdefined.jumps <- TRUE
  if (length(mcmc.jumps$alpha) != 1) illdefined.jumps <- TRUE
  if (length(mcmc.jumps$A) != 1) illdefined.jumps <- TRUE
  if (length(mcmc.jumps$B) != 1) illdefined.jumps <- TRUE
  if (illdefined.jumps) stop("Argument mcmc.jumps is ill-defined!")
  

    
  # Distance and NA matrices
  nsites <- nrow(sites)
  distances <- as.matrix(dist(rbind(sites, knots)))
  dss <- distances[1:nsites, 1:nsites]
  dsk <- distances[1:nsites, -(1:nsites)]
  na.mat <- 1 - is.na(y)  # Matrix of missing values
  
  # Missing values are replaced by 0 but ignored in the MCMC function
  y0 <- y
  y0[is.na(y)] <- 0
  
  
  
  
  ##########################
  ## Calling C++ function ##
  ##########################
  # C++ function
  if (!quiet) cat("MCMC begins...\n")
  
  if (fit.margins) {
    result <- .Call('hkevp_mcmc_hkevp', PACKAGE = 'hkevp', y0, sites, knots, niter, nburn, trace, quiet, dss, dsk, na.mat, spatial.covariates, log.scale, gev.vary, correlation, mcmc.init$loc, mcmc.init$scale, mcmc.init$shape, mcmc.init$range, mcmc.init$sill, mcmc.init$alpha, mcmc.init$tau, mcmc.init$A, mcmc.init$B, mcmc.prior$constant.gev, mcmc.prior$beta.sd, mcmc.prior$range, mcmc.prior$sill, mcmc.prior$alpha, mcmc.prior$tau, mcmc.jumps$gev, mcmc.jumps$range, mcmc.jumps$alpha, mcmc.jumps$tau, mcmc.jumps$A, mcmc.jumps$B)
  }
  else {
    result <- .Call('hkevp_mcmc_deponly', PACKAGE = 'hkevp', y0, sites, knots, niter, nburn, trace, quiet, max(dss), dsk, na.mat, mcmc.init$alpha, mcmc.init$tau, mcmc.init$A, mcmc.init$B, mcmc.prior$alpha, mcmc.prior$tau, mcmc.jumps$alpha, mcmc.jumps$tau, mcmc.jumps$A, mcmc.jumps$B)
  }
  
  
  # Burn-in period and thinning
  mcmc.index <- seq(nburn+1, niter, by = nthin)
  result$GEV <- result$GEV[,,mcmc.index]
  result$alpha <- matrix(result$alpha[mcmc.index])
  result$tau <- matrix(result$tau[mcmc.index])
  result$A <- result$A[,,mcmc.index]
  result$llik <- matrix(result$llik[mcmc.index])
  result$spatial$beta <- result$spatial$beta[mcmc.index,,]
  result$spatial$sills <- result$spatial$sills[mcmc.index,]
  result$spatial$ranges <- result$spatial$ranges[mcmc.index,]
  
  
  if (log.scale) result$GEV[,2,] <- exp(result$GEV[,2,])  # GEV log scale transformed to GEV scale
  result$data <- y
  result$sites <- sites
  result$knots <- knots
  result$spatial.covariates <- spatial.covariates
  result$correlation <- correlation
  result$nstep <- length(mcmc.index)
  result$log.scale <- log.scale
  result$fit.type <- ifelse(fit.margins, "hkevp", "dep-only")
  
  return(result)
}
