#' Conversion of an \code{dmbc_fit} object to an object of class \code{mcmc}.
#' 
#' \code{dmbc_fit_to_mcmc} converts an object of class \code{dmbc_fit}
#'   to an object with class \code{mcmc}.
#' 
#' @param res An object of type \code{dmbc_fit}.
#' @param include.burnin A logical scalar. If \code{TRUE} the burnin
#'   iterations (if available) are not removed.
#' @param verbose A logical scalar. If \code{TRUE} prints additional
#'   warnings during the conversion.
#' @return An object of type \code{mcmc}.
#' @seealso
#'   \code{\link{dmbc}()} for for fitting a DMBC model;
#'   \code{\link{dmbc_fit-class}};
#'   \code{\link[coda]{mcmc}}.
#' @author Sergio Venturini \email{sergio.venturini@unito.it}
#' @examples
#' \dontrun{
#' data(simdiss, package = "dmbc")
#' 
#' G <- 3
#' p <- 2
#' prm.prop <- list(z = 1.5, alpha = .75)
#' burnin <- 2000
#' nsim <- 1000
#' seed <- 2301
#' 
#' set.seed(seed)
#' 
#' control <- list(burnin = burnin, nsim = nsim, z.prop = prm.prop[["z"]],
#'   alpha.prop = prm.prop[["alpha"]], verbose = TRUE)
#' sim.dmbc <- dmbc(simdiss, p, G, control)
#' sim.mcmc <- dmbc_fit_to_mcmc(sim.dmbc@results[[1]], TRUE)
#' plot(sim.mcmc)
#' }
#' @export
dmbc_fit_to_mcmc <- function(res, include.burnin = FALSE, verbose = TRUE) {
  control <- res@control
  burnin <- control[["burnin"]]
  nsim <- control[["nsim"]]
  thin <- control[["thin"]]
  store.burnin <- control[["store.burnin"]]
  totiter <- burnin + nsim

  n <- res@dim[["n"]]
  p <- res@dim[["p"]]
  G <- res@dim[["G"]]
  S <- res@dim[["S"]]

  if (store.burnin) {
    if (include.burnin) {
      tokeep <- seq(1, totiter, by = thin)
      tokeep <- 1:length(tokeep)
    } else {
      todrop <- seq(1, burnin, by = thin)
      tokeep <- seq(1, totiter, by = thin)
      tokeep <- (length(todrop) + 1):length(tokeep)
    }
  } else {
    if (verbose && include.burnin)
      warning("burnin iterations not shown because the 'store.burnin' option was set to FALSE.")
    tokeep <- seq(1, nsim, by = thin)
    tokeep <- 1:length(tokeep)
  }

  theta <- matrix(NA, nrow = length(tokeep), ncol = (2*n*p*G + 4*G + 2*S*G + S + 3))
  theta_nm <- character(2*n*p*G + 4*G + 2*S*G + S + 3)
  count_z <- count_g <- 1
  for (g in 1:G) {
    for (j in 1:p) {
      for (i in 1:n) {
        theta[, count_z] <- res@z.chain[tokeep, i, j, g]
        theta_nm[count_z] <- paste0("z[", i, ", ", j, ", ", g, "]")
        theta[, n*p*G + count_z] <- res@z.chain.p[tokeep, i, j, g]
        theta_nm[n*p*G + count_z] <- paste0("z_p[", i, ", ", j, ", ", g, "]")
        count_z <- count_z + 1
      }
    }
    theta[, 2*n*p*G + count_g] <- res@alpha.chain[tokeep, g]
    theta_nm[2*n*p*G + count_g] <- paste0("alpha[", g, "]")
    theta[, 2*n*p*G + G + count_g] <- res@eta.chain[tokeep, g]
    theta_nm[2*n*p*G + G + count_g] <- paste0("eta[", g, "]")
    theta[, 2*n*p*G + 2*G + count_g] <- res@sigma2.chain[tokeep, g]
    theta_nm[2*n*p*G + 2*G + count_g] <- paste0("sigma2[", g, "]")
    theta[, 2*n*p*G + 3*G + count_g] <- res@lambda.chain[tokeep, g]
    theta_nm[2*n*p*G + 3*G + count_g] <- paste0("lambda[", g, "]")
    for (s in 1:S) {
      theta[, 2*n*p*G + 4*G + S*(g - 1) + s] <- res@prob.chain[tokeep, s, g]
      theta_nm[2*n*p*G + 4*G + S*(g - 1) + s] <- paste0("prob[", s, g, "]")
      theta[, 2*n*p*G + 4*G + S*G + S*(g - 1) + s] <- res@x.ind.chain[tokeep, s, g]
      theta_nm[2*n*p*G + 4*G + S*G + S*(g - 1) + s] <- paste0("x_ind[", s, g, "]")
    }
    count_g <- count_g + 1
  }
  for (s in 1:S) {
    theta[, 2*n*p*G + 4*G + 2*S*G + s] <- res@x.chain[tokeep, s]
    theta_nm[2*n*p*G + 4*G + 2*S*G + s] <- paste0("x[", s, "]")
  }
  theta[, 2*n*p*G + 4*G + 2*S*G + S + 1] <- res@dens$loglik[tokeep]
  theta_nm[2*n*p*G + 4*G + 2*S*G + S + 1] <- paste0("loglik")
  theta[, 2*n*p*G + 4*G + 2*S*G + S + 2] <- res@dens$logprior[tokeep]
  theta_nm[2*n*p*G + 4*G + 2*S*G + S + 2] <- paste0("logprior")
  theta[, 2*n*p*G + 4*G + 2*S*G + S + 3] <- res@dens$logpost[tokeep]
  theta_nm[2*n*p*G + 4*G + 2*S*G + S + 3] <- paste0("logpost")
  
  colnames(theta) <- theta_nm

  if (store.burnin) {
    if (include.burnin) {
      out <- coda::mcmc(theta, start = 1, end = totiter, thin = thin)
    } else {
      out <- coda::mcmc(theta, start = (burnin + 1), end = totiter, thin = thin)
    }
  } else {
    out <- coda::mcmc(theta, start = (burnin + 1), end = totiter, thin = thin)
  }

  return(out)
}

#' Conversion of an \code{dmbc_fit_list} object to a \code{list}.
#' 
#' \code{dmbc_fit_list_to_list} converts an object of class
#'   \code{dmbc_fit_list} to a list of arrays including all the parameter.
#'   chains. It is intended for internal use mainly.
#' 
#' @param res An object of type \code{dmbc_fit_list}.
#' @param include.burnin A logical scalar. If \code{TRUE} the burnin
#'   iterations (if available) are not removed.
#' @param verbose A logical scalar. If \code{TRUE} prints additional
#'   warnings during the conversion.
#' @return An object of type \code{mcmc.list}.
#' @seealso
#'   \code{\link{dmbc}()} for for fitting a DMBC model;
#'   \code{\link{dmbc_fit_list-class}}.
#' @author Sergio Venturini \email{sergio.venturini@unito.it}
#' @examples
#' \dontrun{
#' data(simdiss, package = "dmbc")
#' 
#' G <- 3
#' p <- 2
#' prm.prop <- list(z = 1.5, alpha = .75)
#' burnin <- 2000
#' nsim <- 1000
#' seed <- 2301
#' 
#' set.seed(seed)
#' 
#' control <- list(burnin = burnin, nsim = nsim, z.prop = prm.prop[["z"]],
#'   alpha.prop = prm.prop[["alpha"]], nchains = 2, verbose = TRUE)
#' sim.dmbc <- dmbc(simdiss, p, G, control)
#' sim.list <- dmbc_fit_list_to_list(sim.dmbc, TRUE)
#' 
#' library(bayesplot)
#' mcmc_trace(sim.list, regex_pars = "lambda")
#' }
#' @export
dmbc_fit_list_to_list <- function(res, include.burnin = FALSE, verbose = TRUE) {
  control <- res@results[[1]]@control
  burnin <- control[["burnin"]]
  nsim <- control[["nsim"]]
  thin <- control[["thin"]]
  nchains <- control[["nchains"]]
  store.burnin <- control[["store.burnin"]]
  totiter <- burnin + nsim

  n <- res@results[[1]]@dim[["n"]]
  p <- res@results[[1]]@dim[["p"]]
  G <- res@results[[1]]@dim[["G"]]
  S <- res@results[[1]]@dim[["S"]]

  if (store.burnin) {
    if (include.burnin) {
      tokeep <- seq(1, totiter, by = thin)
      tokeep <- 1:length(tokeep)
    } else {
      todrop <- seq(1, burnin, by = thin)
      tokeep <- seq(1, totiter, by = thin)
      tokeep <- (length(todrop) + 1):length(tokeep)
    }
  } else {
    if (verbose && include.burnin)
      warning("burnin iterations not shown because the 'store.burnin' option was set to FALSE.")
    tokeep <- seq(1, nsim, by = thin)
    tokeep <- 1:length(tokeep)
  }

  out <- list()
  for (c in 1:nchains) {
    theta <- matrix(NA, nrow = length(tokeep), ncol = (2*n*p*G + 4*G + 2*S*G + S + 3))
    theta_nm <- character(2*n*p*G + 4*G + 2*S*G + S + 3)
    count_z <- count_g <- 1
    for (g in 1:G) {
      for (j in 1:p) {
        for (i in 1:n) {
          theta[, count_z] <- res@results[[c]]@z.chain[tokeep, i, j, g]
          theta_nm[count_z] <- paste0("z[", i, ", ", j, ", ", g, "]")
          theta[, n*p*G + count_z] <- res@results[[c]]@z.chain.p[tokeep, i, j, g]
          theta_nm[n*p*G + count_z] <- paste0("z_p[", i, ", ", j, ", ", g, "]")
          count_z <- count_z + 1
        }
      }
      theta[, 2*n*p*G + count_g] <- res@results[[c]]@alpha.chain[tokeep, g]
      theta_nm[2*n*p*G + count_g] <- paste0("alpha[", g, "]")
      theta[, 2*n*p*G + G + count_g] <- res@results[[c]]@eta.chain[tokeep, g]
      theta_nm[2*n*p*G + G + count_g] <- paste0("eta[", g, "]")
      theta[, 2*n*p*G + 2*G + count_g] <- res@results[[c]]@sigma2.chain[tokeep, g]
      theta_nm[2*n*p*G + 2*G + count_g] <- paste0("sigma2[", g, "]")
      theta[, 2*n*p*G + 3*G + count_g] <- res@results[[c]]@lambda.chain[tokeep, g]
      theta_nm[2*n*p*G + 3*G + count_g] <- paste0("lambda[", g, "]")
      for (s in 1:S) {
        theta[, 2*n*p*G + 4*G + S*(g - 1) + s] <- res@results[[c]]@prob.chain[tokeep, s, g]
        theta_nm[2*n*p*G + 4*G + S*(g - 1) + s] <- paste0("prob[", s, g, "]")
        theta[, 2*n*p*G + 4*G + S*G + S*(g - 1) + s] <- res@results[[c]]@x.ind.chain[tokeep, s, g]
        theta_nm[2*n*p*G + 4*G + S*G + S*(g - 1) + s] <- paste0("x_ind[", s, g, "]")
      }
      count_g <- count_g + 1
    }
    for (s in 1:S) {
      theta[, 2*n*p*G + 4*G + 2*S*G + s] <- res@results[[c]]@x.chain[tokeep, s]
      theta_nm[2*n*p*G + 4*G + 2*S*G + s] <- paste0("x[", s, "]")
    }
    theta[, 2*n*p*G + 4*G + 2*S*G + S + 1] <- res@results[[c]]@dens$loglik[tokeep]
    theta_nm[2*n*p*G + 4*G + 2*S*G + S + 1] <- paste0("loglik")
    theta[, 2*n*p*G + 4*G + 2*S*G + S + 2] <- res@results[[c]]@dens$logprior[tokeep]
    theta_nm[2*n*p*G + 4*G + 2*S*G + S + 2] <- paste0("logprior")
    theta[, 2*n*p*G + 4*G + 2*S*G + S + 3] <- res@results[[c]]@dens$logpost[tokeep]
    theta_nm[2*n*p*G + 4*G + 2*S*G + S + 3] <- paste0("logpost")
    
    colnames(theta) <- theta_nm
    out[[c]] <- theta
  }

  return(out)
}

#' Conversion of an \code{dmbc_fit_list} object to an object of class
#'   \code{mcmc.list}.
#' 
#' \code{dmbc_fit_list_to_mcmc.list} converts an object of class
#'   \code{dmbc_fit_list} to an object with class \code{mcmc.list}.
#' 
#' @param res An object of type \code{dmbc_fit_list}.
#' @param include.burnin A logical scalar. If \code{TRUE} the burnin
#'   iterations (if available) are not removed.
#' @param verbose A logical scalar. If \code{TRUE} prints additional
#'   warnings during the conversion.
#' @return An object of type \code{mcmc.list}.
#' @seealso
#'   \code{\link{dmbc}()} for for fitting a DMBC model;
#'   \code{\link{dmbc_fit_list-class}};
#'   \code{\link[coda]{mcmc.list}}.
#' @author Sergio Venturini \email{sergio.venturini@unito.it}
#' @examples
#' \dontrun{
#' data(simdiss, package = "dmbc")
#' 
#' G <- 3
#' p <- 2
#' prm.prop <- list(z = 1.5, alpha = .75)
#' burnin <- 2000
#' nsim <- 1000
#' seed <- 2301
#' 
#' set.seed(seed)
#' 
#' control <- list(burnin = burnin, nsim = nsim, z.prop = prm.prop[["z"]],
#'   alpha.prop = prm.prop[["alpha"]], nchains = 2, verbose = TRUE)
#' sim.dmbc <- dmbc(simdiss, p, G, control)
#' sim.mcmc <- dmbc_fit_list_to_mcmc.list(sim.dmbc, TRUE)
#' plot(sim.mcmc)
#' }
#' @export
dmbc_fit_list_to_mcmc.list <- function(res, include.burnin = FALSE, verbose = TRUE) {
  control <- res@results[[1]]@control
  burnin <- control[["burnin"]]
  nsim <- control[["nsim"]]
  thin <- control[["thin"]]
  nchains <- control[["nchains"]]
  store.burnin <- control[["store.burnin"]]
  totiter <- burnin + nsim

  out <- dmbc_fit_list_to_list(res, include.burnin = include.burnin, verbose = verbose)
  if (store.burnin) {
    if (include.burnin) {
      out <- lapply(out, coda::mcmc, start = 1, end = totiter, thin = thin)
    } else {
      out <- lapply(out, coda::mcmc, start = (burnin + 1), end = totiter, thin = thin)
    }
  } else {
    out <- lapply(out, coda::mcmc, start = (burnin + 1), end = totiter, thin = thin)
  }
  out <- coda::mcmc.list(out)

  return(out)
}
