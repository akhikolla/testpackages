#' Fitter function for DMBC models.
#'
#' \code{dmbc_fit()} is the main function that estimates a DMBC model.
#'
#' @param D A list whose elements are the dissimilarity matrices corresponding
#'   to the judgments expressed by the \emph{S} subjects/raters. These matrices
#'   must be defined as a \code{dist} object.
#' @param p A length-one numeric vector indicating the number of dimensions of the
#'   latent space.
#' @param G A length-one numeric vector indicating the number of cluster to
#'   partition the \emph{S} subjects.
#' @param family A length-one character vector representing the type of data to
#'   analyze. Currently, it accepts only the 'binomial' value, but future
#'   developments will include the possibility to analyze continuous,
#'   multinomial and count data.
#' @param control A list of control parameters that affect the sampling
#'   but do not affect the posterior distribution See
#'   \code{\link{dmbc_control}()} for more details.
#' @param prior A list containing the prior hyperparameters. See
#'   \code{\link{dmbc_prior}()} for more details.
#' @param start A named list of starting values for the MCMC algorithm (see
#'   \code{\link{dmbc_init}}).
#' @return A \code{dmbc_fit_list} object.
#' @author Sergio Venturini \email{sergio.venturini@unito.it}
#' @seealso \code{\link{dmbc_data}} for a description of the data format.
#' @seealso \code{\link{dmbc_fit_list}} for a description of the elements
#'   included in the returned object.
#' @references
#'   Venturini, S., Piccarreta, R. (2019), "A Bayesian Approach for Model-Based
#'   Clustering of Several Binary Dissimilarity Matrices: the \pkg{dmbc}
#'   Package in \code{R}", Technical report.
#' @examples
#' \dontrun{
#' data(simdiss, package = "dmbc")
#'
#' G <- 3
#' p <- 2
#' prm.prop <- list(z = 1.5, alpha = .75)
#' burnin <- 20000
#' nsim <- 10000
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
#' summary(sim.dmbc, include.burnin = FALSE)
#'
#' library(bayesplot)
#' library(ggplot2)
#' color_scheme_set("teal")
#' plot(sim.dmbc, what = "trace", regex_pars = "eta")
#'
#' z <- dmbc_get_configuration(sim.dmbc, chain = 1, est = "mean",
#'   labels = 1:16)
#' summary(z)
#' color_scheme_set("mix-pink-blue")
#' graph <- plot(z, size = 2, size_lbl = 3)
#' graph + panel_bg(fill = "gray90", color = NA)
#' }
#' @export
dmbc_fit <- function(D, p, G, family, control, prior, start) {
	S <- length(D)
	n <- attr(D[[1]], "Size")
	m <- n*(n - 1)/2
	totiter <- control[["burnin"]] + control[["nsim"]]
	p <- as.integer(p)
	G <- as.integer(G)
	
	z.chain <- z.chain.p <- array(NA, dim = c(totiter, n, p, G))
	eta.chain <- alpha.chain <- sigma2.chain <- lambda.chain <- array(NA, dim = c(totiter, G))
	x.chain <- array(NA, dim = c(totiter, S))
	prob.chain <- x.ind.chain <- array(0, dim = c(totiter, S, G))
	loglik <- logprior <- logpost <- numeric(totiter)
	
	# recover prior hyperparameters
  hyper.eta.a <- prior[["eta"]][["a"]]
  hyper.eta.b <- prior[["eta"]][["b"]]
	hyper.sigma2.a <- prior[["sigma2"]][["a"]]
	hyper.sigma2.b <- prior[["sigma2"]][["b"]]
	hyper.lambda <- prior[["lambda"]]
	
	# start iteration
	if (control[["verbose"]]) message("Running the MCMC simulation...")
	
	res.mcmc <- .Call('dmbc_mcmc', PACKAGE = 'dmbc',
		raiD = as.integer(unlist(D)),
		raix = as.integer(start$x),
		raing = as.integer(start$ng),
		radalpha = as.double(start$alpha),
		rn = as.integer(n),
		rp = as.integer(p),
		rG = as.integer(G),
		rS = as.integer(S),
		rtotiter = as.integer(totiter),
		radZ = as.double(start$z),
		rgamma_z = as.double(control[["z.prop"]]),
		reta = as.double(start$eta),
		rgamma_alpha = as.double(control[["alpha.prop"]]),
		rsigma2 = as.double(start$sigma2),
		rlambda = as.double(start$lambda),
		rhyper_eta_a = as.double(hyper.eta.a),
		rhyper_eta_b = as.double(hyper.eta.b),
    rhyper_sigma2_a = as.double(hyper.sigma2.a),
    rhyper_sigma2_b = as.double(hyper.sigma2.b),
		rhyper_lambda = as.double(hyper.lambda),
    rfamily = as.character(family),
		rverbose = as.integer(control[["verbose"]])
	)

	z.chain <- z.chain.p <- array(res.mcmc[[1]], c(totiter, n, p, G))
	alpha.chain <- array(res.mcmc[[2]], c(totiter, G))
	eta.chain <- array(res.mcmc[[3]], c(totiter, G))
	sigma2.chain <- array(res.mcmc[[4]], c(totiter, G))
	lambda.chain <- array(res.mcmc[[5]], c(totiter, G))
	prob.chain <- array(res.mcmc[[6]], c(totiter, S, G))
	x.chain <- array(res.mcmc[[7]], c(totiter, S))
	x.ind.chain <- array(res.mcmc[[8]], c(totiter, S, G))
	accept <- t(array(res.mcmc[[9]], c(G, 2)))
	loglik <- as.numeric(res.mcmc[[10]])
	logprior <- as.numeric(res.mcmc[[11]])
	logpost <- as.numeric(res.mcmc[[12]])

  if (control[["procrustes"]] | control[["relabel"]]) {
  	# post-processing:
  	if (control[["verbose"]]) message("Post-processing the chain:")

    if (control[["procrustes"]]) {
    	## Procrustes transformation of Z_g
    	if (control[["verbose"]]) message("   - applying Procrustes transformation...")
      if (control[["verbose"]]) {
        pb <- dmbc_pb(min = 0, max = (totiter*G - 1), width = 49)
      }
      no <- 0
    	for (niter in 1:totiter) {
    		for (g in 1:G) {
          if (control[["verbose"]]) dmbc_setpb(pb, no)
    			if (p == 1) {
    				z.chain.p[niter, , , g] <- as.numeric(MCMCpack::procrustes(as.matrix(z.chain[niter, , , g]),
              as.matrix(z.chain[totiter, , , g]), translation = TRUE, dilation = FALSE)$X.new)
    			} else {
    				z.chain.p[niter, , , g] <- MCMCpack::procrustes(z.chain[niter, , , g], z.chain[totiter, , , g],
              translation = TRUE, dilation = FALSE)$X.new
    			}
          no <- no + 1
    		}
    	}
      if (control[["verbose"]]) {
        # message("done!")
        close(pb)
      }
    }

    if (control[["relabel"]]) {
    	# relabel the parameter chain
    	if (G > 1) {
    		if (totiter > 10) {
          if (control[["verbose"]]) message("   - relabeling the parameter chain...")
    			init <- ifelse(totiter <= 100, 5, 100)
    			
    			theta <- .Call('dmbc_pack_par', PACKAGE = 'dmbc',
    				radz = as.double(z.chain.p),
    				radalpha = as.double(alpha.chain),
    				radlambda = as.double(lambda.chain),
    				rn = as.integer(n),
    				rp = as.integer(p),
    				rM = as.integer(totiter),
    				rG = as.integer(G)
    			)

    			theta.relab <- .Call('dmbc_relabel', PACKAGE = 'dmbc',
    				radtheta = as.double(theta),
    				radz = as.double(z.chain.p),
    				radalpha = as.double(alpha.chain),
    				radeta = as.double(eta.chain),
    				radsigma2 = as.double(sigma2.chain),
    				radlambda = as.double(lambda.chain),
    				radprob = as.double(prob.chain),
    				raix_ind = as.integer(x.ind.chain),
    				rinit = as.integer(init),
    				rn = as.integer(n),
    				rp = as.integer(p),
    				rS = as.integer(S),
    				rM = as.integer(totiter),
    				rR = as.integer(m + 1),
    				rG = as.integer(G),
            rverbose = as.integer(control[["verbose"]])
    			)

    			theta <- array(theta.relab[[1]], c(totiter, (m + 1), G))  # this is not needed elsewhere
    			z.chain.p <- array(theta.relab[[2]], c(totiter, n, p, G))
    			alpha.chain <- array(theta.relab[[3]], c(totiter, G))
    			eta.chain <- array(theta.relab[[4]], c(totiter, G))
    			sigma2.chain <- array(theta.relab[[5]], c(totiter, G))
    			lambda.chain <- array(theta.relab[[6]], c(totiter, G))
    			prob.chain <- array(theta.relab[[7]], c(totiter, S, G))
    			x.ind.chain <- array(theta.relab[[8]], c(totiter, S, G))
    			x.chain <- t(apply(x.ind.chain, 1, function(x) as.integer(x %*% 1:G)))

    			# if (control[["verbose"]]) # message("done!")
    		} else {
    			warning("the number of iterations is too small for relabeling; relabeling skipped.", call. = FALSE,
            immediate. = TRUE)
    		}
    	}
    }
  }

  # apply thinning
  if (control[["thin"]] > 1) {
    if (control[["store.burnin"]]) {
      tokeep <- seq(1, totiter, by = control[["thin"]])
    } else {
      tokeep <- seq(control[["burnin"]] + 1, totiter, by = control[["thin"]])
    }
  } else {
    if (control[["store.burnin"]]) {
      tokeep <- seq(1, totiter, by = 1)
    } else {
      tokeep <- seq(control[["burnin"]] + 1, totiter, by = 1)
    }
  }
  z.chain <- z.chain[tokeep, , , , drop = FALSE]
  z.chain.p <- z.chain.p[tokeep, , , , drop = FALSE]
  alpha.chain <- alpha.chain[tokeep, , drop = FALSE]
  eta.chain <- eta.chain[tokeep, , drop = FALSE]
  sigma2.chain <- sigma2.chain[tokeep, , drop = FALSE]
  lambda.chain <- lambda.chain[tokeep, , drop = FALSE]
  prob.chain <- prob.chain[tokeep, , , drop = FALSE]
  x.ind.chain <- x.ind.chain[tokeep, , , drop = FALSE]
  x.chain <- x.chain[tokeep, , drop = FALSE]
  loglik <- loglik[tokeep]
  logprior <- logprior[tokeep]
  logpost <- logpost[tokeep]

  # return results
	out <- new("dmbc_fit",
		z.chain = z.chain,
		z.chain.p = z.chain.p,
		alpha.chain = alpha.chain,
		eta.chain = eta.chain,
		sigma2.chain = sigma2.chain,
		lambda.chain = lambda.chain,
		prob.chain = prob.chain,
		x.ind.chain = x.ind.chain,
		x.chain = x.chain,
		accept = accept,
		diss = D,
		dens = list(loglik = loglik, logprior = logprior, logpost = logpost),
    control = control,
    prior = prior,
		dim = list(n = n, p = p, G = G, S = S),
    model = new("dmbc_model", p = p, G = G, family = family)
	)

	return(out)
}

#' Log-likelihood for DMBC models.
#'
#' \code{dmbc_logLik_rbmds()} computes the log-likelihood value for a DMBC model.
#'
#' @param D A list whose elements are the dissimilarity matrices corresponding
#'   to the judgments expressed by the \emph{S} subjects/raters. These matrices
#'   must be defined as a \code{dist} object.
#' @param Z A numeric matrix containing the latent configuration.
#' @param alpha A numeric vector containing the alpha values.
#'
#' @return A length-one numeric vector of the log-likelihood value.
#'
#' @author Sergio Venturini \email{sergio.venturini@unito.it}
#'
#' @seealso \code{\link{dmbc}()}.
#'
#' @references
#'   Venturini, S., Piccarreta, R. (2019), "A Bayesian Approach for Model-Based
#'   Clustering of Several Binary Dissimilarity Matrices: the \pkg{dmbc}
#'   Package in \code{R}", Technical report.
#'
#' @export
dmbc_logLik_rbmds <- function(D, Z, alpha) {
	S <- length(D)
	diss <- list.sum(D)
	delta <- dist(Z)
	pi <- expit(alpha + delta)
	# out <- sum(diss*log(pi) + (S - diss)*log(1 - pi))
	out <- sum(log((pi^diss)*((1 - pi)^(S - diss))))
	
	return(out)
}

#' Log-likelihood for DMBC models.
#'
#' \code{dmbc_logLik()} computes the log-likelihood value for a DMBC model.
#'
#' @param D A list whose elements are the dissimilarity matrices corresponding
#'   to the judgments expressed by the \emph{S} subjects/raters. These matrices
#'   must be defined as a \code{dist} object.
#' @param Z A numeric matrix containing the latent configuration.
#' @param alpha A numeric vector containing the alpha values.
#' @param lambda A numeric vector containing the alpha lambda.
#' @param x A numeric vector containing the cluster indicator values.
#'
#' @return A length-one numeric vector of the log-likelihood value.
#'
#' @author Sergio Venturini \email{sergio.venturini@unito.it}
#'
#' @seealso \code{\link{dmbc}()}.
#'
#' @references
#'   Venturini, S., Piccarreta, R. (2019), "A Bayesian Approach for Model-Based
#'   Clustering of Several Binary Dissimilarity Matrices: the \pkg{dmbc}
#'   Package in \code{R}", Technical report.
#'
#' @export
dmbc_logLik <- function(D, Z, alpha, lambda, x) {
	G <- length(alpha)
	ng <- as.numeric(table(factor(x, levels = 1:G)))
	logfg <- numeric(G)
	ll <- 0
	for (g in 1:G) {
		if (ng[g] > 0) {
			xg <- which(x == g)
			logfg[g] <- dmbc_logLik_rbmds(D[xg], Z[, , g], alpha[g])
			ll <- ll + ng[g]*log(lambda[g]) + logfg[g]
		}
	}

	return(ll)
}
