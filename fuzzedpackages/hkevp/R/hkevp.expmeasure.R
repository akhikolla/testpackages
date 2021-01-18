#' @title
#' Exponent measure of the HKEVP
#' 
#' @author
#' Quentin Sebille
#' 
#' 
#' @description 
#' Exponent measure \eqn{V(z_1,...,z_n)} of the HKEVP of \cite{Reich and Shaby (2012)}, with given model parameters or output from \code{hkevp.fit} or \code{latent.fit}.
#' 
#' 
#' @param z
#' The vector \eqn{(z_1,...,z_n)} where the exponent measure is computed. Can be of length one and thus corresponds then to \eqn{(z,...,z)}.
#' 
#' @param fit
#' Output from the \code{hkevp.fit} procedure.
#' 
#' @param sites
#' The coordinates of the sites where the data are observed. Each row correspond to a site position.
#' 
#' @param knots
#' The coordinates of the knots in the HKEVP. By default, the positions of the knots coincide with the positions of the sites.
#' 
#' @param alpha
#' The dependence parameter \eqn{\alpha} of the HKEVP: a single value in (0,1].
#' 
#' @param tau
#' The bandwidth parameter \eqn{\tau} of the kernel functions in the HKEVP: a positive value.
#' 
#' 
#' 
#' 
#' 
#' @details
#' The exponent measure describes the spatial dependence structure of a max-stable process, independently from the values of the marginal parameters. If \eqn{Z(\cdot)} is a simple max-stable process, i.e. with unit GEV(1,1,1) margins, recorded at the set of sites \eqn{(s_1, \ldots, s_n)}, its joint cumulative probability density function is given by:
#' \deqn{P\{ Z(s_1)\leq z_1, \ldots, Z(s_n)\leq z_n \} = \exp(-V(z_1, \ldots, z_n)) ~,}
#' where \eqn{V} is the so-called exponent measure.

#' For the HKEVP, the exponent measure is explicit for any number \eqn{n} of sites:
#' \deqn{V(z_1, \ldots, z_n) = \sum_{\ell=1}^L \left[ \sum_{i=1}^n \left(\frac{\omega_\ell(s_i)}{z_i}\right)^{1/\alpha}\right]^{\alpha} ~.}
#' 
#' If argument \code{fit} is provided, the predictive distribution of \deqn{V(z_1, \ldots, z_n)} is computed. If not, the function uses arguments \code{sites}, \code{knots}, \code{alpha}, and \code{tau}.
#' 
#' 
#'
#' @return
#' Either a vector if argument \code{fit} is provided, or a single value.
#' 
#' 
#' 
#' @export
#' 
#' 
#' @references 
#' Reich, B. J., & Shaby, B. A. (2012). A hierarchical max-stable spatial model for extreme precipitation. The annals of applied statistics, 6(4), 1430. <DOI:10.1214/12-AOAS591>
#' 
#'
#' 
#' 
#' 
#' @examples
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
#' predict.em <- hkevp.expmeasure(1, fit = fit)
#' true.em <- hkevp.expmeasure(1, sites, sites, alpha, tau)
#' # plot(predict.em, ylim = range(predict.em, true.em), type = "l")
#' # abline(h = true.em, col = 2, lwd = 2)
#' 
#' 
hkevp.expmeasure <- function(z, sites, knots, alpha, tau, fit) {
  if (missing(fit)) {  ## Evaluation with given parameters
    # Catch Errors and default values
    if (alpha <= 0 | alpha > 1) stop("alpha must be between 0 and 1!")
    if (tau <= 0) stop("tau must be positive!")
    if (length(z) == 1) z <- rep(z, nrow(sites))
    if (length(z) != nrow(sites)) stop("z and sites do not match!")
    
    # Computing the kernel matrix omega
    nsites <- nrow(sites)
    if (nsites == 1) return(1/z)
    dsk <- as.matrix(dist(rbind(sites, knots)))[1:nsites, -(1:nsites)]
    omega <- exp(-dsk^2/(2*tau^2))
    omega <- sweep(omega, MARGIN = 1, STATS = rowSums(omega), FUN = "/")
    
    # Computing the exponent measure
    result <- sum(colSums(sweep(x = omega, MARGIN = 1, STATS = z, FUN = '/') ^ (1/alpha) ) ^ alpha)
    return(result)
  } else {
    sites <- fit$sites
    knots <- fit$knots
    alpha <- fit$alpha
    tau <- fit$tau
    
    if (length(z) == 1) z <- rep(z, nrow(sites))
    if (length(z) != nrow(sites)) stop("Arguments z and sites do not match!")
    nsites <- nrow(sites)
    if (nsites == 1) return(1/z)
    
    if (fit$fit.type == "hkevp" | fit$fit.type == "dep-only") {
      dsk <- as.matrix(dist(rbind(sites, knots)))[1:nsites, -(1:nsites)]
      fun1 <- function(param, z, dsk) {
        omega <- exp(-dsk^2/(2*param[2]^2))
        omega <- sweep(omega, MARGIN = 1, STATS = rowSums(omega), FUN = "/")
        sum(colSums(sweep(x = omega, MARGIN = 1, STATS = z, FUN = '/') ^ (1/param[1]) ) ^ param[1])
      }
      return(apply(cbind(alpha, tau), 1, fun1, z, dsk))
    }
    else if (fit$fit.type == "latent") {
      result <- rep(sum(1/z), fit$nstep)
    }
    else stop("Incorrect fit.type!")
  }
}



