#' Model selection of DMBC models.
#'
#' \code{dmbc_IC()} is the main function for simultaneously selecting the
#'   optimal latent space dimension (\emph{p}) and number of clusters
#'   (\emph{G}) for a DMBC analysis.
#'
#' @param data An object of class \code{dmbc_data} containing the data
#'   to analyze.
#' @param pmax A length-one numeric vector indicating the maximum number of
#'   dimensions of the latent space to consider.
#' @param Gmax A length-one numeric vector indicating the maximum number of
#'   cluster to consider.
#' @param control A list of control parameters that affect the sampling
#'   but do not affect the posterior distribution See
#'   \code{\link{dmbc_control}()} for more details.
#' @param prior A list containing the prior hyperparameters. See
#'   \code{\link{dmbc_prior}()} for more details.
#' @param est A length-one character vector indicating the estimate type to
#'   use. Possible values are \code{mean}, \code{median}, \code{ml} and
#'   \code{map}.
#'
#' @return A \code{dmbc_ic} object.
#'
#' @author Sergio Venturini \email{sergio.venturini@unito.it}
#'
#' @seealso \code{\link{dmbc}()} for fitting a DMBC model.
#' @seealso \code{\link{dmbc_ic}} for a description of the elements included
#'   in the returned object.
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
#' pmax <- 2
#' Gmax <- 2
#' prm.prop <- list(z = 1.5, alpha = .75)
#' burnin <- 2000
#' nsim <- 1000
#' seed <- 1809
#' 
#' set.seed(seed)
#' 
#' control <- list(burnin = burnin, nsim = nsim, z.prop = prm.prop[["z"]],
#'   alpha.prop = prm.prop[["alpha"]], random.start = TRUE, verbose = TRUE,
#'   thin = 10, store.burnin = TRUE)
#' sim.ic <- dmbc_IC(data = simdiss, pmax = pmax, Gmax = Gmax, control = control,
#'   est = "mean")
#' 
#' pmax <- pmax + 1
#' Gmax <- Gmax + 2
#' new.ic <- update(sim.ic, pmax = pmax, Gmax = Gmax)
#' new.ic
#' 
#' # plot the results
#' library(bayesplot)
#' library(ggplot2)
#' color_scheme_set("mix-yellow-blue")
#' p <- plot(new.ic, size = c(4, 1.5))
#' p + panel_bg(fill = "gray90", color = NA)
#' }
#' @export
dmbc_IC <- function(data, pmax = 3, Gmax = 5, control = dmbc_control(), prior = NULL, est = "mean") {
  D <- data@diss
  control <- check_list_na(control, dmbc_control())
  if (!check_control(control))
    stop("the control list is not correct; see the documentation for more details.")
  if (control[["nchains"]] > 1)
    stop("the number of chains to simulate must be set to 1 when computing the DMBC information criterion.")
  if (!(est %in% c("mean", "median", "ml", "map"))) {
    stop("the estimate type specified is not available.")
  }
  
  verbose <- control[["verbose"]]
  logprior <- logmlik <- logcorrfact <- dcic <- matrix(NA, nrow = pmax, ncol = Gmax)
	res_list <- list()
	res_save <- list()
	res_all <- list()

	res.i <- 1
	for (G.i in 1:Gmax) {
		for (p.i in 1:pmax) {
			if (verbose)
			  message("--- p = ", p.i, " -- G = ", G.i, " ---")
			
			res <- dmbc(data = data, p = p.i, G = G.i, control = control, prior = prior)
			res_list[[p.i]] <- res
			if (p.i == pmax)
				res_save[[G.i]] <- res
			
      if (est == "mean") {
        est.tmp <- dmbc_get_postmean(res)
			} else if (est == "median") {
				est.tmp <- dmbc_get_postmedian(res)
			} else if (est == "ml") {
				est.tmp <- dmbc_get_ml(res)
			} else if (est == "map") {
				est.tmp <- dmbc_get_map(res)
			}
      z.m <- est.tmp$z
      alpha.m <- est.tmp$alpha
      eta.m <- est.tmp$eta
      sigma2.m <- est.tmp$sigma2
      lambda.m <- est.tmp$lambda
			class.m <- dmbc_get_postmean(res)$cluster
			res_all[[res.i]] <- list(z.m = z.m, alpha.m = alpha.m, eta.m = eta.m, sigma2.m = sigma2.m, lambda.m = lambda.m, class.m = class.m)
			names(res_all)[res.i] <- paste("p = ", p.i, " -- G = ", G.i, sep = "")

			res.i <- res.i + 1
			
			logprior[p.i, G.i] <- log_marg_prior(res, z.m)
			logmlik[p.i, G.i] <- log_marg_lik(res, z.m)
			if (p.i > 1)
				logcorrfact[p.i, G.i] <- log_corr_fact(res_list[[p.i - 1]], z.m)
			dcic[p.i, G.i] <-  -2*(logprior[p.i, G.i] + logmlik[p.i, G.i] +
        ifelse(p.i > 1, sum(logcorrfact[2:p.i, G.i], na.rm = TRUE), 0))
			
			if (verbose) {
				message(" ")
			}
		}
		res_list <- list()
	}
	
	out <- new("dmbc_ic",
		logprior = logprior,
		logmlik = logmlik,
		logcorrfact = logcorrfact,
		DCIC = dcic,
		post.est = res_all,
		est = est,
		res_last_p = res_save
	)

	return(out)
}

log_marg_lik <- function(res, Z) {
	D <- res@results[[1]]@diss
	control <- res@results[[1]]@control
  prior <- res@results[[1]]@prior

  burnin <- control[["burnin"]]
  nsim <- control[["nsim"]]
  thin <- control[["thin"]]
  store.burnin <- control[["store.burnin"]]
  totiter <- burnin + nsim

  if (store.burnin) {
    todrop <- seq(1, burnin, by = thin)
    tokeep <- seq(1, totiter, by = thin)
    tokeep <- (length(todrop) + 1):length(tokeep)
  } else {
    tokeep <- seq(1, nsim, by = thin)
    tokeep <- 1:length(tokeep)
  }

	x <- dmbc_get_postmean(res)$cluster
	n <- dim(Z)[1]
	p <- dim(Z)[2]
	G <- dim(Z)[3]
	if (G > 1) {
		q <- 3*G
		theta <- cbind(res@results[[1]]@alpha.chain[tokeep, ], res@results[[1]]@sigma2.chain[tokeep, ],
      res@results[[1]]@lambda.chain[tokeep, ])
		theta.star <- pcaPP::l1median(theta)
		theta.star[(2*G + 1):(3*G)] <- colMeans(theta)[(2*G + 1):(3*G)]   # needed cause the 'ddirichlet' function
                                                                      # returns a 0 when the sum of lambdas != 1
		H.star <- robustbase::covMcd(theta[, -q])$cov   # last dimension removed otherwise the hessian is singular
	} else {
		q <- 2
		theta <- cbind(res@results[[1]]@alpha.chain[tokeep, 1], res@results[[1]]@sigma2.chain[tokeep, 1])
		theta.star <- pcaPP::l1median(theta)
		H.star <- robustbase::covMcd(theta)$cov
	}
	if (G > 1) {
		loglik <- dmbc_logLik(D, Z, theta.star[1:G], theta.star[(2*G + 1):(3*G)], x)
	} else {
		loglik <- dmbc_logLik(D, Z, theta.star[1], theta.star[2], x)
	}
	alpha <- prior[["sigma2"]][["a"]]
	beta <- prior[["sigma2"]][["b"]]
	if (G > 1) {
		lambda.hyp <- prior[["lambda"]]
	}
	logprior <- 0
	for (g in 1:G) {
		logprior <- logprior + dnorm(theta.star[g], sd = sqrt(theta.star[G + g]), log = TRUE)
		logprior <- logprior + dinvgamma(theta.star[G + g], alpha = alpha, beta = beta, log = TRUE)
	}
	if (G > 1) {
		logprior <- logprior + log(ddirichlet(theta.star[(2*G + 1):(3*G)], lambda.hyp))
	}
	logmlik <- q*log(2*pi)/2 + log(det(H.star))/2 + loglik + logprior

	return(logmlik)
}

log_marg_prior <- function(res, Z) {
	n <- dim(Z)[1]
	p <- dim(Z)[2]
	G <- dim(Z)[3]
	a_g <- res@results[[1]]@prior[["eta"]][["a"]]
	b_g <- res@results[[1]]@prior[["eta"]][["b"]]
	logprior <- 0
	for (g in 1:G) {
		logprior <- logprior + lgamma(a_g[g] + n*p/2) - lgamma(a_g[g]) + a_g[g]*log(b_g[g]) - (a_g[g] +
      n*p/2)*log(b_g[g] + sum(Z[, , g]^2)/2)
	}
	logprior <- logprior - n*p*G*log(2*pi)/2

	return(logprior)
}

log_corr_fact <- function(res, Z) {
	n <- dim(Z)[1]
	p <- dim(Z)[2]
	a <- res@results[[1]]@prior[["eta"]][["a"]][1]
	b <- res@results[[1]]@prior[["eta"]][["b"]][1]
	logcorrfact <- n*log(2*pi)/2 + lgamma(a + n*p/2) - lgamma(a + n*(p + 1)/2) + n*log(b + sum(Z^2)/2)/2

	return(logcorrfact)
}
