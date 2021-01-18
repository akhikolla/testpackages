#' Extractor function for a fitted DMBC model.
#'
#' \code{dmbc_get_postmean()} is an extractor function for extracting the
#'   posterior mean estimates of the parameters for a fitted DMBC model.
#'
#' @param res An object of class \code{dmbc_fit_list}.
#' @param chain A length-one numeric vector indicating the MCMC chain number
#'   to use.
#'
#' @return A named \code{list} with the following elements:
#'   \describe{
#'     \item{\code{z}: }{array of latent coordinates posterior mean estimates}
#'     \item{\code{alpha}: }{numeric vector of alpha posterior mean estimates}
#'     \item{\code{eta}: }{numeric vector of eta posterior mean estimates}
#'     \item{\code{sigma2}: }{numeric vector of sigma2 posterior mean estimates}
#'     \item{\code{lambda}: }{numeric vector of lambda posterior mean estimates}
#'     \item{\code{prob}: }{numeric matrix of probability posterior mean estimates}
#'     \item{\code{cluster}: }{numeric vector of cluster membership posterior
#'       mean estimates}
#'     \item{\code{chain}: }{length-one numeric vector of the MCMC chain number
#'       used}
#'   }
#'
#' @author Sergio Venturini \email{sergio.venturini@unito.it}
#'
#' @seealso \code{\link{dmbc_data}} for a description of the data format.
#' @seealso \code{\link{dmbc_fit_list}} for a description of the elements
#'   included in the returned object.
#'
#' @references
#'   Venturini, S., Piccarreta, R. (2019), "A Bayesian Approach for Model-Based
#'   Clustering of Several Binary Dissimilarity Matrices: the \pkg{dmbc}
#'   Package in \code{R}", Technical report.
#'
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
#'   alpha.prop = prm.prop[["alpha"]], random.start = TRUE, verbose = TRUE,
#'   nchains = 2, thin = 10, store.burnin = TRUE, threads = 2,
#'   parallel = "snow")
#' sim.dmbc <- dmbc(simdiss, p, G, control)
#'
#' dmbc_get_postmean(sim.dmbc, chain = 1)
#' }
#' @export
dmbc_get_postmean <- function(res, chain = 1) {
  res_chain <- res@results[[chain]]
  control <- res@results[[1]]@control
  burnin <- control[["burnin"]]
  nsim <- control[["nsim"]]
  thin <- control[["thin"]]
  nchains <- control[["nchains"]]
  store.burnin <- control[["store.burnin"]]
  totiter <- burnin + nsim
  
  if (chain > nchains)
    stop("the specified chain is not available.")

  if (store.burnin) {
    todrop <- seq(1, burnin, by = thin)
    tokeep <- seq(1, totiter, by = thin)
    tokeep <- (length(todrop) + 1):length(tokeep)
  } else {
    tokeep <- seq(1, nsim, by = thin)
    tokeep <- 1:length(tokeep)
  }

  z.postmean <- apply(res_chain@z.chain.p[tokeep, , , , drop = FALSE], c(2, 3, 4), mean, na.rm = TRUE)
  alpha.postmean <- colMeans(res_chain@alpha.chain[tokeep, , drop = FALSE], na.rm = TRUE)
  eta.postmean <- colMeans(res_chain@eta.chain[tokeep, , drop = FALSE], na.rm = TRUE)
  sigma2.postmean <- colMeans(res_chain@sigma2.chain[tokeep, , drop = FALSE], na.rm = TRUE)
  lambda.postmean <- colMeans(res_chain@lambda.chain[tokeep, , drop = FALSE], na.rm = TRUE)
  prob.postmean <- apply(res_chain@prob.chain[tokeep, , , drop = FALSE], c(2, 3), mean, na.rm = TRUE)
  class.postmean <- apply(prob.postmean, 1, which.max)
  out <- list(z = z.postmean,
              alpha = alpha.postmean,
              eta = eta.postmean,
              sigma2 = sigma2.postmean,
              lambda = lambda.postmean,
              prob = prob.postmean,
              cluster = class.postmean,
              chain = chain)

  return(out)
}

#' Extractor function for a fitted DMBC model.
#'
#' \code{dmbc_get_postmedian()} is an extractor function for extracting the
#'   posterior median estimates of the parameters for a fitted DMBC model.
#'
#' @param res An object of class \code{dmbc_fit_list}.
#' @param chain A length-one numeric vector indicating the MCMC chain number
#'   to use.
#'
#' @return A named \code{list} with the following elements:
#'   \describe{
#'     \item{\code{z}: }{array of latent coordinates posterior median estimates}
#'     \item{\code{alpha}: }{numeric vector of alpha posterior median estimates}
#'     \item{\code{eta}: }{numeric vector of eta posterior median estimates}
#'     \item{\code{sigma2}: }{numeric vector of sigma2 posterior median estimates}
#'     \item{\code{lambda}: }{numeric vector of lambda posterior median estimates}
#'     \item{\code{prob}: }{numeric matrix of probability posterior median estimates}
#'     \item{\code{cluster}: }{numeric vector of cluster membership posterior
#'       median estimates}
#'     \item{\code{chain}: }{length-one numeric vector of the MCMC chain number
#'       used}
#'   }
#'
#' @author Sergio Venturini \email{sergio.venturini@unito.it}
#'
#' @seealso \code{\link{dmbc_data}} for a description of the data format.
#' @seealso \code{\link{dmbc_fit_list}} for a description of the elements
#'   included in the returned object.
#'
#' @references
#'   Venturini, S., Piccarreta, R. (2019), "A Bayesian Approach for Model-Based
#'   Clustering of Several Binary Dissimilarity Matrices: the \pkg{dmbc}
#'   Package in \code{R}", Technical report.
#'
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
#'   alpha.prop = prm.prop[["alpha"]], random.start = TRUE, verbose = TRUE,
#'   nchains = 2, thin = 10, store.burnin = TRUE, threads = 2,
#'   parallel = "snow")
#' sim.dmbc <- dmbc(simdiss, p, G, control)
#'
#' dmbc_get_postmedian(sim.dmbc, chain = 1)
#' }
#' @export
dmbc_get_postmedian <- function(res, chain = 1) {
  res_chain <- res@results[[chain]]
  control <- res@results[[1]]@control
  burnin <- control[["burnin"]]
	nsim <- control[["nsim"]]
  thin <- control[["thin"]]
  nchains <- control[["nchains"]]
  store.burnin <- control[["store.burnin"]]
	totiter <- burnin + nsim
	
  if (chain > nchains)
    stop("the specified chain is not available.")

  if (store.burnin) {
    todrop <- seq(1, burnin, by = thin)
    tokeep <- seq(1, totiter, by = thin)
    tokeep <- (length(todrop) + 1):length(tokeep)
  } else {
    tokeep <- seq(1, nsim, by = thin)
    tokeep <- 1:length(tokeep)
  }

	z.postmedian <- apply(res_chain@z.chain.p[tokeep, , , , drop = FALSE], c(2, 3, 4), median, na.rm = TRUE)
	alpha.postmedian <- colMedians(res_chain@alpha.chain[tokeep, , drop = FALSE], na.rm = TRUE)
	eta.postmedian <- colMedians(res_chain@eta.chain[tokeep, , drop = FALSE], na.rm = TRUE)
	sigma2.postmedian <- colMedians(res_chain@sigma2.chain[tokeep, , drop = FALSE], na.rm = TRUE)
	lambda.postmedian <- colMedians(res_chain@lambda.chain[tokeep, , drop = FALSE], na.rm = TRUE)
	prob.postmedian <- apply(res_chain@prob.chain[tokeep, , , drop = FALSE], c(2, 3), median, na.rm = TRUE)
	class.postmedian <- apply(prob.postmedian, 1, which.max)
	out <- list(z = z.postmedian,
              alpha = alpha.postmedian,
              eta = eta.postmedian,
              sigma2 = sigma2.postmedian,
		          lambda = lambda.postmedian,
              prob = prob.postmedian,
              cluster = class.postmedian,
              chain = chain)

	return(out)
}

#' Extractor function for a fitted DMBC model.
#'
#' \code{dmbc_get_ml()} is an extractor function for extracting the
#'   maximum likelihood estimates of the parameters for a fitted DMBC model.
#'
#' @param res An object of class \code{dmbc_fit_list}.
#' @param chain A length-one numeric vector indicating the MCMC chain number
#'   to use.
#'
#' @return A named \code{list} with the following elements:
#'   \describe{
#'     \item{\code{z}: }{array of latent coordinates posterior mean estimates}
#'     \item{\code{alpha}: }{numeric vector of alpha posterior mean estimates}
#'     \item{\code{eta}: }{numeric vector of eta posterior mean estimates}
#'     \item{\code{sigma2}: }{numeric vector of sigma2 posterior mean estimates}
#'     \item{\code{lambda}: }{numeric vector of lambda posterior mean estimates}
#'     \item{\code{prob}: }{numeric matrix of probability posterior mean estimates}
#'     \item{\code{cluster}: }{numeric vector of cluster membership posterior
#'       mean estimates}
#'     \item{\code{loglik}: }{length-one numeric vector of the maximum
#'       log-likelihood value}
#'     \item{\code{chain}: }{length-one numeric vector of the MCMC chain number
#'       used}
#'   }
#'
#' @author Sergio Venturini \email{sergio.venturini@unito.it}
#'
#' @seealso \code{\link{dmbc_data}} for a description of the data format.
#' @seealso \code{\link{dmbc_fit_list}} for a description of the elements
#'   included in the returned object.
#'
#' @references
#'   Venturini, S., Piccarreta, R. (2019), "A Bayesian Approach for Model-Based
#'   Clustering of Several Binary Dissimilarity Matrices: the \pkg{dmbc}
#'   Package in \code{R}", Technical report.
#'
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
#'   alpha.prop = prm.prop[["alpha"]], random.start = TRUE, verbose = TRUE,
#'   nchains = 2, thin = 10, store.burnin = TRUE, threads = 2,
#'   parallel = "snow")
#' sim.dmbc <- dmbc(simdiss, p, G, control)
#'
#' dmbc_get_ml(sim.dmbc, chain = 1)
#' }
#' @export
dmbc_get_ml <- function(res, chain = 1) {
  res_chain <- res@results[[chain]]
  control <- res@results[[1]]@control
  burnin <- control[["burnin"]]
  nsim <- control[["nsim"]]
  thin <- control[["thin"]]
  nchains <- control[["nchains"]]
  store.burnin <- control[["store.burnin"]]
  totiter <- burnin + nsim
	ll <- res_chain@dens$loglik
  
  if (chain > nchains)
    stop("the specified chain is not available.")

  if (store.burnin) {
    todrop <- seq(1, burnin, by = thin)
    tokeep <- seq(1, totiter, by = thin)
    tokeep <- (length(todrop) + 1):length(tokeep)
  } else {
    tokeep <- seq(1, nsim, by = thin)
    tokeep <- 1:length(tokeep)
  }

  i.ml <- which.max(ll[tokeep])
  if (store.burnin) {
    i.ml <- length(todrop) + i.ml
  }

  prob.ml <- res_chain@prob.chain[i.ml, , ]
	out <- list(z = res_chain@z.chain.p[i.ml, , , ],
              alpha = res_chain@alpha.chain[i.ml, ],
          		eta = res_chain@eta.chain[i.ml, ],
              sigma2 = res_chain@sigma2.chain[i.ml, ],
              lambda = res_chain@lambda.chain[i.ml, ],
              prob = prob.ml,
          		cluster = apply(prob.ml, 1, which.max),
              loglik = ll[i.ml],
              chain = chain)

	return(out)
}

#' Extractor function for a fitted DMBC model.
#'
#' \code{dmbc_get_map()} is an extractor function for extracting the
#'   maximum-a-posterior estimates of the parameters for a fitted DMBC model.
#'
#' @param res An object of class \code{dmbc_fit_list}.
#' @param chain A length-one numeric vector indicating the MCMC chain number
#'   to use.
#'
#' @return A named \code{list} with the following elements:
#'   \describe{
#'     \item{\code{z}: }{array of latent coordinates posterior mean estimates}
#'     \item{\code{alpha}: }{numeric vector of alpha posterior mean estimates}
#'     \item{\code{eta}: }{numeric vector of eta posterior mean estimates}
#'     \item{\code{sigma2}: }{numeric vector of sigma2 posterior mean estimates}
#'     \item{\code{lambda}: }{numeric vector of lambda posterior mean estimates}
#'     \item{\code{prob}: }{numeric matrix of probability posterior mean estimates}
#'     \item{\code{cluster}: }{numeric vector of cluster membership posterior
#'       mean estimates}
#'     \item{\code{logpost}: }{length-one numeric vector of the maximum
#'       log-posterior value}
#'     \item{\code{chain}: }{length-one numeric vector of the MCMC chain number
#'       used}
#'   }
#'
#' @author Sergio Venturini \email{sergio.venturini@unito.it}
#'
#' @seealso \code{\link{dmbc_data}} for a description of the data format.
#' @seealso \code{\link{dmbc_fit_list}} for a description of the elements
#'   included in the returned object.
#'
#' @references
#'   Venturini, S., Piccarreta, R. (2019), "A Bayesian Approach for Model-Based
#'   Clustering of Several Binary Dissimilarity Matrices: the \pkg{dmbc}
#'   Package in \code{R}", Technical report.
#'
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
#'   alpha.prop = prm.prop[["alpha"]], random.start = TRUE, verbose = TRUE,
#'   nchains = 2, thin = 10, store.burnin = TRUE, threads = 2,
#'   parallel = "snow")
#' sim.dmbc <- dmbc(simdiss, p, G, control)
#'
#' dmbc_get_map(sim.dmbc, chain = 1)
#' }
#' @export
dmbc_get_map <- function(res, chain = 1) {
  res_chain <- res@results[[chain]]
  control <- res@results[[1]]@control
  burnin <- control[["burnin"]]
  nsim <- control[["nsim"]]
  thin <- control[["thin"]]
  store.burnin <- control[["store.burnin"]]
  nchains <- control[["nchains"]]
  totiter <- burnin + nsim
	lpost <- res_chain@dens$logpost
  
  if (chain > nchains)
    stop("the specified chain is not available.")

  if (store.burnin) {
    todrop <- seq(1, burnin, by = thin)
    tokeep <- seq(1, totiter, by = thin)
    tokeep <- (length(todrop) + 1):length(tokeep)
  } else {
    tokeep <- seq(1, nsim, by = thin)
    tokeep <- 1:length(tokeep)
  }

  i.map <- which.max(lpost[tokeep])
  if (store.burnin) {
    i.map <- length(todrop) + i.map
  }

	prob.map <- res_chain@prob.chain[i.map, , ]
  out <- list(z = res_chain@z.chain.p[i.map, , , ],
              alpha = res_chain@alpha.chain[i.map, ],
          		eta = res_chain@eta.chain[i.map, ],
              sigma2 = res_chain@sigma2.chain[i.map, ],
          		lambda = res_chain@lambda.chain[i.map, ],
              prob = prob.map,
              cluster = apply(prob.map, 1, which.max),
              logpost = lpost[i.map],
              chain = chain)

	return(out)
}

#' Extractor function for a fitted DMBC model.
#'
#' \code{dmbc_get_configuration()} is an extractor function for extracting the
#'   latent configuration estimates of a fitted DMBC model.
#'
#' @param res An object of class \code{dmbc_fit_list}.
#' @param chain A length-one numeric vector indicating the MCMC chain number
#'   to use.
#' @param est A length-one character vector indicating the estimate type to use.
#' @param labels An optional character vector with the object labels.
#'
#' @return A \code{\link{dmbc_config}} object.
#'
#' @author Sergio Venturini \email{sergio.venturini@unito.it}
#'
#' @seealso \code{\link{dmbc_data}} for a description of the data format.
#' @seealso \code{\link{dmbc_fit_list}} for a description of the elements
#'   included in the returned object.
#'
#' @references
#'   Venturini, S., Piccarreta, R. (2019), "A Bayesian Approach for Model-Based
#'   Clustering of Several Binary Dissimilarity Matrices: the \pkg{dmbc}
#'   Package in \code{R}", Technical report.
#'
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
#'   alpha.prop = prm.prop[["alpha"]], random.start = TRUE, verbose = TRUE,
#'   nchains = 2, thin = 10, store.burnin = TRUE, threads = 2,
#'   parallel = "snow")
#' sim.dmbc <- dmbc(simdiss, p, G, control)
#'
#' z <- dmbc_get_configuration(sim.dmbc, chain = 1, est = "mean")
#' summary(z)
#'
#' library(bayesplot)
#' library(ggplot2)
#' color_scheme_set("mix-pink-blue")
#' graph <- plot(z, size = 2, size_lbl = 3)
#' graph + panel_bg(fill = "gray90", color = NA)
#' }
#' @export
dmbc_get_configuration <- function(res, chain = 1, est = "mean", labels = character(0)) {
  res_chain <- res@results[[chain]]
  control <- res@results[[1]]@control
  burnin <- control[["burnin"]]
  nsim <- control[["nsim"]]
  thin <- control[["thin"]]
  nchains <- control[["nchains"]]
  store.burnin <- control[["store.burnin"]]
  totiter <- burnin + nsim
  
  if (chain > nchains)
    stop("the specified chain is not available.")
  if (!(est %in% c("mean", "median", "ml", "map")))
    stop("the estimate type specified is not available.")
  labels <- as.character(labels)
  if (length(labels) && (length(labels) != res_chain@dim[["n"]]))
    stop("the number of labels provided must be equal to the number of objects in the data.")

  if (store.burnin) {
    todrop <- seq(1, burnin, by = thin)
    tokeep <- seq(1, totiter, by = thin)
    tokeep <- (length(todrop) + 1):length(tokeep)
  } else {
    tokeep <- seq(1, nsim, by = thin)
    tokeep <- 1:length(tokeep)
  }

  res.est <- switch(est,
    mean = dmbc_get_postmean(res, chain = chain),
    median = dmbc_get_postmedian(res, chain = chain),
    ml = dmbc_get_ml(res, chain = chain),
    map = dmbc_get_map(res, chain = chain))
  Z.est <- res.est$z
  Z.sd <- apply(res_chain@z.chain.p[tokeep, , , , drop = FALSE], c(2, 3, 4), sd, na.rm = TRUE)
  dimnames(Z.est)[[1]] <- dimnames(Z.sd)[[1]] <- if (length(labels)) labels else paste0("i = ", 1:dim(Z.est)[1])
  dimnames(Z.est)[[2]] <- dimnames(Z.sd)[[2]] <- paste0("p_", 1:dim(Z.est)[2])
  dimnames(Z.est)[[3]] <- dimnames(Z.sd)[[3]] <- paste0("g = ", 1:dim(Z.est)[3])
  cl <- res.est$cluster

  out <- new("dmbc_config",
    Z.est = Z.est,
    Z.sd = Z.sd,
    cluster = cl,
    est = est,
    n = res_chain@dim[["n"]],
    p = res_chain@dim[["p"]],
    S = res_chain@dim[["S"]],
    G = res_chain@dim[["G"]],
    family = res_chain@model@family,
    chain = chain,
    labels = labels)

  return(out)
}

#' Auxiliary function for checking the grouping results of a fitted DMBC model.
#'
#' \code{dmbc_check_groups()} is an auxiliary function for checking whether
#'   the cluster membership estimates provided by the individual chains of the
#'   fitted model provided agree or not.
#'
#' @param res An object of class \code{dmbc_fit_list}.
#' @param est A length-one character vector indicating the estimate type to use.
#'
#' @return A length-one logical vector, which is equal to TRUE if all simulated chains
#'   provide the same cluster membership estimates, and FALSE otherwise.
#'
#' @author Sergio Venturini \email{sergio.venturini@unito.it}
#'
#' @seealso \code{\link{dmbc_get_configuration}()} for a description of the
#'   configuration extractor function.
#' @seealso \code{\link{dmbc_fit_list}} for a description of a fitted
#'   DMBC model.
#'
#' @references
#'   Venturini, S., Piccarreta, R. (2019), "A Bayesian Approach for Model-Based
#'   Clustering of Several Binary Dissimilarity Matrices: the \pkg{dmbc}
#'   Package in \code{R}", Technical report.
#'
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
#'   alpha.prop = prm.prop[["alpha"]], random.start = TRUE, verbose = TRUE,
#'   nchains = 2, thin = 10, store.burnin = TRUE, threads = 2,
#'   parallel = "snow")
#' sim.dmbc <- dmbc(simdiss, p, G, control)
#'
#' dmbc_check_groups(sim.dmbc)
#' }
#' 
#' @importFrom stats chisq.test
#' @export
dmbc_check_groups <- function(res, est = "mean") {
  consistent <- TRUE
  control <- res@results[[1]]@control
  dims <- res@results[[1]]@dim
  nchains <- control$nchains
  S <- dims[["S"]]
  G <- dims[["G"]]
  cluster <- matrix(NA, nrow = S, ncol = nchains)
  cluster_tbl <- array(NA, dim = c(G, G, nchains*(nchains - 1)/2))
  cluster_chk <- logical(nchains*(nchains - 1)/2)
  cluster_count <- 1

  suppressWarnings(
    if (nchains > 1) {
      for (ch in 1:nchains) {
        cluster[, ch] <- clusters(dmbc_get_configuration(res, est = est))
      }
      for (i in 1:(nchains - 1)) {
        for (j in (i + 1):nchains) {
          cluster_tbl[, , cluster_count] <- table(factor(cluster[, i], levels = 1:G),
                                                  factor(cluster[, j], levels = 1:G))
          cluster_chk[cluster_count] <- all.equal(stats::chisq.test(cluster_tbl[, , cluster_count])$statistic,
            S*(G - 1), check.attributes = FALSE)
          cluster_count <- cluster_count + 1
        }
      }
      if (!all(cluster_chk)) consistent <- FALSE
    }
  )
  
  return(consistent)
}

#' Auxiliary function for realigning the grouping of a fitted DMBC model.
#'
#' \code{dmbc_match_groups()} is an auxiliary function for realigning the
#'   cluster membership estimates provided by the individual chains of the
#'   fitted model if they do not agree.
#'
#' @param res An object of class \code{dmbc_fit_list}.
#' @param est A length-one character vector indicating the estimate type to use.
#' @param ref A length-one numeric vector indicating the chain number to use as
#'   the reference.
#'
#' @return An object of class \code{dmbc_fit_list}.
#'
#' @author Sergio Venturini \email{sergio.venturini@unito.it}
#'
#' @seealso \code{\link{dmbc_check_groups}()} for checking the consistency
#'   of the cluster memberships across chains for a fitted DMBC model.
#' @seealso \code{\link{dmbc_get_configuration}()} for a description of the
#'   configuration extractor function.
#' @seealso \code{\link{dmbc_fit_list}} for a description of a fitted
#'   DMBC model.
#'
#' @references
#'   Venturini, S., Piccarreta, R. (2019), "A Bayesian Approach for Model-Based
#'   Clustering of Several Binary Dissimilarity Matrices: the \pkg{dmbc}
#'   Package in \code{R}", Technical report.
#'
#' @examples
#' \dontrun{
#' data(simdiss, package = "dmbc")
#'
#' G <- 5
#' p <- 3
#' prm.prop <- list(z = 4, alpha = 2)
#' burnin <- 2000
#' nsim <- 1000
#' seed <- 2301
#'
#' set.seed(seed)
#'
#' control <- list(burnin = burnin, nsim = nsim, z.prop = prm.prop[["z"]],
#'   alpha.prop = prm.prop[["alpha"]], random.start = TRUE, verbose = TRUE,
#'   nchains = 6, store.burnin = TRUE, threads = 2, parallel = "snow")
#' sim.dmbc <- dmbc(simdiss, p, G, control)
#'
#' sim.dmbc_new <- dmbc_match_groups(sim.dmbc)
#' }
#' 
#' @importFrom stats chisq.test
#' @export
dmbc_match_groups <- function(res, est = "mean", ref = 1) {
  control <- res@results[[1]]@control
  dims <- res@results[[1]]@dim
  nchains <- control$nchains
  S <- dims[["S"]]
  G <- dims[["G"]]
  cluster <- matrix(NA, nrow = S, ncol = nchains)
  cluster_tbl <- array(NA, dim = c(G, G, nchains))

  if (nchains > 1) {
    for (ch in 1:nchains) {
      cluster[, ch] <- clusters(dmbc_get_configuration(res, est = est, chain = ch))
    }
    for (i in (1:nchains)[-ref]) {
      cluster_tbl[, , i] <- table(factor(cluster[, i], levels = 1:G),
                                  factor(cluster[, ref], levels = 1:G))
      chisq <- stats::chisq.test(cluster_tbl[, , i])$statistic
      if (is.nan(chisq)) {
        # [[TODO: the following condition is not general enough]]
        cluster_chk <- all(sort(margin.table(cluster_tbl[, , i], 1)) == 
          sort(margin.table(cluster_tbl[, , i], 2)))
      } else {
        cluster_chk <- all.equal(chisq, S*(G - 1), check.attributes = FALSE)
      }
      if (!cluster_chk) {
        warning("the cluster memberships of some chains do not fully agree.")
        return(res)
      } else {
        if (sum(diag(cluster_tbl[, , i])) != S) {
          new_cluster <- apply(cluster_tbl[, , i], 2, which.max)
          res@results[[i]]@z.chain <- res@results[[i]]@z.chain[, , , new_cluster, drop = FALSE]
          res@results[[i]]@z.chain.p <- res@results[[i]]@z.chain.p[, , , new_cluster, drop = FALSE]
          res@results[[i]]@alpha.chain <- res@results[[i]]@alpha.chain[, new_cluster, drop = FALSE]
          res@results[[i]]@eta.chain <- res@results[[i]]@eta.chain[, new_cluster, drop = FALSE]
          res@results[[i]]@sigma2.chain <- res@results[[i]]@sigma2.chain[, new_cluster, drop = FALSE]
          res@results[[i]]@lambda.chain <- res@results[[i]]@lambda.chain[, new_cluster, drop = FALSE]
          res@results[[i]]@prob.chain <- res@results[[i]]@prob.chain[, , new_cluster, drop = FALSE]
          for (it in 1:dim(res@results[[i]]@x.chain)[1]) {
            x <- x_tmp <- res@results[[i]]@x.chain[it, ]
            for (g in 1:G) {
              x <- replace(x, x_tmp == g, new_cluster[g])
            }
            res@results[[i]]@x.chain[it, ] <- x
          }
          res@results[[i]]@x.ind.chain <- res@results[[i]]@x.ind.chain[, , new_cluster, drop = FALSE]
          res@results[[i]]@accept <- res@results[[i]]@accept[, new_cluster, drop = FALSE]
        }
      }
    }
  }

  return(res)
}
