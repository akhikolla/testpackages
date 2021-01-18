#' Sum of squared residuals (SSR) from the observed distances and the given
#' latent configuration.
#' 
#' \code{comp_ssr} computes the sum of squared residuals (SSR) from the
#'   observed distances (\code{diss}) and the given latent coordinates
#'   (\code{x}).
#' 
#' @param x Real matrix containing the latent configuration.
#' @param diss Observed dissimilarities (provided as a distance matrix).
#' @return A length-one numeric vector providing the SSR for its arguments.
#' @seealso \code{\link{bmds}} for (one-way) Bayesian (metric) multidimensional
#'   scaling.
#' @author Sergio Venturini \email{sergio.venturini@unito.it}
#' @examples
#' n <- 10000
#' nr <- 200
#' nc <- floor(n/nr)
#' x <- matrix(rnorm(1:n), nrow = nr, ncol = nc)
#' obsdiss <- dist(x)
#' ssr <- numeric(ncol(x))
#' for (i in 1:ncol(x)) {
#'   ssr[i] <- comp_ssr(x[, 1:i], obsdiss)
#' }
#' plot(ssr, xlab = "number of dimensions", ylab = "SSR", type = "b")
#' @export
comp_ssr <- function(x, diss) {
	d <- as.matrix(diss)
	delta <- as.matrix(dist(x))
	ssr <- sum((d[lower.tri(d)] - delta[lower.tri(delta)])^2)
	
	return(ssr)
}

#' Adjustment of the center and orientation of a latent configuration.
#' 
#' \code{adjust_x} adjusts the center and orientation of a latent configuration
#'   in Bayesian (metric) multidimensional scaling (BMDS).
#' 
#' @param x Numeric matrix containing the latent configuration.
#' @seealso \code{\link{bmds}} for (one-way) Bayesian (metric) multidimensional
#'   scaling.
#' @return A list with elements:
#' \describe{
#'   \item{\code{x}}{A real matrix containing the adjusted latent
#'                   configuration.}
#'   \item{\code{Sig_x}}{The variance and covariance matrix of the adjusted
#'                       latent configuration.}
#' }
#' @author Sergio Venturini \email{sergio.venturini@unito.it}
#' @examples
#' n <- 100
#' nr <- 20
#' nc <- floor(n/nr)
#' x <- matrix(rnorm(1:n), nrow = nr, ncol = nc)
#' adj_x <- adjust_x(x)
#' adj_x$x
#' adj_x$Sig_x
#' @export
adjust_x <- function(x) {
	x_new <- scale(x, center = TRUE, scale = FALSE)
	n <- dim(x)[1]
	p <- dim(x)[2]
	Sig_x <- (t(x_new) %*% x_new)/n
	
	e <- eigen(Sig_x, symmetric = TRUE)
	eval <- e$values[seq_len(p)]
	evec <- e$vectors[, seq_len(p), drop = FALSE]
	rotat <- t(evec)
	
	x_new <- t(rotat %*% t(x_new))
	
	Sig_x <- (rotat %*% Sig_x) %*% t(rotat)
	
	return(list(x = x_new, Sig_x = Sig_x))
}

#' Posterior mode latent configuration in Bayesian multidimensional scaling (BMDS).
#' 
#' \code{bmds_get_x_mode} returns the latent configuration that produced the
#'   largest posterior value during the MCMC.
#' 
#' @param D Observed dissimilarities (provided as a distance matrix).
#' @param res Results of a BMDS analysis as obtained with the
#'   \code{\link{bmds}} function.
#' @param p.i A length-one numeric vector providing the index of the solution to
#'   use.
#' @param min_p A length-one numeric vector providing the minimum value of the
#'   latent space dimension to use.
#' @param max_p A length-one numeric vector providing the maximum value
#'   of the latent space dimension to use.
#' @param start A length-one numeric vector providing the iteration
#'   number to start from.
#' @param end A length-one numeric vector providing the iteration
#'   number where to end.
#' @return A real matrix containing the posterior mode latent configuration.
#' @seealso \code{\link{bmds}} for (one-way) Bayesian (metric) multidimensional
#'   scaling.
#' @author Sergio Venturini \email{sergio.venturini@unito.it}
#' @examples
#' \dontrun{
#' # Airline Distances Between Cities
#' airline <- read.csv(file = system.file("extdata", "airline.csv",
#'   package = "dmbc"))
#' airline.nm <- airline[, 1]
#' airline <- airline[, 2:31]
#' colnames(airline) <- airline.nm
#' airline <- as.dist(airline)
#' 
#' min_p <- 1
#' max_p <- 4
#' burnin <- 200
#' nsim <- 1000
#' totiter <- burnin + nsim
#' 
#' airline.mds <- cmdscale(airline, max_p)
#' airline.bmds <- bmds(airline, min_p, max_p, burnin, nsim)
#' 
#' opar <- par(mfrow = c(1, 2))
#' plot(min_p:max_p, airline.bmds$mdsIC$mdsic, type = "b",
#'   main = "MDS Information Criterion", xlab = "p", ylab = "MDSIC")
#' MDSICmin <- which.min(airline.bmds$mdsIC$mdsic)
#' points((min_p:max_p)[MDSICmin], airline.bmds$mdsIC$mdsic[MDSICmin],
#'   col = "red", pch = 10, cex = 1.75, lwd = 1.5)
#' 
#' airline.bmds.x.mode <- bmds_get_x_mode(airline, airline.bmds, MDSICmin,
#'   min_p, max_p, start = (burnin + 1), end = totiter)
#' airline.bmds.d <- dist(airline.bmds.x.mode)
#' airline.mds.d <- dist(airline.mds[, 1:((min_p:max_p)[MDSICmin])])
#' plot(airline, airline.bmds.d, type = "n", xlab = "observed",
#'   ylab = "estimated", main = "Airline Distances \n Between Cities",
#'   xlim = c(0, max(airline, airline.bmds.d)),
#'   ylim = c(0, max(airline, airline.bmds.d)))
#' abline(0, 1, lty = 2, col = "gray")
#' points(airline, airline.mds.d, pch = 19, col = "cyan", cex = .5)
#' points(airline, airline.bmds.d, pch = 19, col = "magenta", cex = .5)
#' legend(x = "bottomright", legend = c("Classical MDS", "Bayesian MDS"),
#'   pch = c(19, 19), col = c("cyan", "magenta"))
#' par(opar)
#' }
#' @export
bmds_get_x_mode <- function(D, res, p.i, min_p, max_p, start, end) {
	bmds.x <- res$x.chain[start:end, p.i, , 1:((min_p:max_p)[p.i]), drop = FALSE]
	bmds.ssr <- drop(apply(bmds.x, c(1, 2), comp_ssr, diss = D))
	bmds.x.mode <- drop(bmds.x[which.min(bmds.ssr), , , ])
	
	return(bmds.x.mode)
}

#' Bayesian multidimensional scaling (BMDS) using Markov Chain Monte Carlo
#'   (MCMC).
#' 
#' \code{bmds} computes the Bayesian multidimensional scaling (BMDS) solutions
#'   using Markov Chain Monte Carlo for a range of specified latent space
#'   dimensions.
#' 
#' @param D Observed dissimilarities (provided as a distance matrix).
#' @param min_p A length-one numeric vector providing the minimum value
#'   of the latent space dimension to use.
#' @param max_pm1 A length-one numeric vector providing the maximum
#'   value of the latent space dimension to use (minus 1).
#' @param burnin A length-one numeric vector providing the number of
#'   iterations to use for burnin.
#' @param nsim A length-one numeric vector providing the number of
#'   iterations to use in the MCMC simulation after burnin.
#' @param ic Logical scalar. If \code{TRUE} computes the MDS
#'   information criterion (MDSIC) for all solution requested.
#' @param verbose Logical scalar. If \code{TRUE} prints information
#'   regarding the evolution of the simulation.
#' @return A list with the following elements:
#'   \describe{
#'     \item{\code{x.chain}}{MCMC chain of the latent configuration
#'       coordinates.}
#'     \item{\code{sigma.chain}}{MCMC chain of the random error.}
#'     \item{\code{lambda.chain}}{MCMC chain of the latent configuration
#'       variances.}
#'     \item{\code{stress}}{Numeric vector of the stress function values.}
#'     \item{\code{mdsIC}}{List with two elements, the MDSIC and BIC values
#'       for the required solutions.}
#'     \item{\code{accept}}{Numeric matrix of acceptance rates.}
#'   }
#' @seealso \code{\link{cmdscale}} for classical (metric) multidimensional scaling.
#' @references
#'   Oh, M.-S., Raftery, A. E. (2001), "Bayesian Multidimensional Scaling and
#'   Choice of Dimension", Journal of the American Statistical Association,
#'   96, 1031-1044.
#' @author Sergio Venturini \email{sergio.venturini@unito.it}
#' @examples
#' \dontrun{
#' # Airline Distances Between Cities
#' airline <- read.csv(file = system.file("extdata", "airline.csv",
#'   package = "dmbc"))
#' airline.nm <- airline[, 1]
#' airline <- airline[, 2:31]
#' colnames(airline) <- airline.nm
#' airline <- as.dist(airline)
#' 
#' min_p <- 1
#' max_p <- 4
#' burnin <- 200
#' nsim <- 1000
#' totiter <- burnin + nsim
#' 
#' airline.mds <- cmdscale(airline, max_p)
#' airline.bmds <- bmds(airline, min_p, max_p, burnin, nsim)
#' 
#' opar <- par(mfrow = c(1, 2))
#' plot(min_p:max_p, airline.bmds$mdsIC$mdsic, type = "b",
#'   main = "MDS Information Criterion", xlab = "p", ylab = "MDSIC")
#' MDSICmin <- which.min(airline.bmds$mdsIC$mdsic)
#' points((min_p:max_p)[MDSICmin], airline.bmds$mdsIC$mdsic[MDSICmin],
#'   col = "red", pch = 10, cex = 1.75, lwd = 1.5)
#' 
#' airline.bmds.x.mode <- bmds_get_x_mode(airline, airline.bmds, MDSICmin,
#'   min_p, max_p, start = (burnin + 1), end = totiter)
#' airline.bmds.d <- dist(airline.bmds.x.mode)
#' airline.mds.d <- dist(airline.mds[, 1:((min_p:max_p)[MDSICmin])])
#' plot(airline, airline.bmds.d, type = "n", xlab = "observed",
#'   ylab = "estimated", main = "Airline Distances \n Between Cities",
#'   xlim = c(0, max(airline, airline.bmds.d)),
#'   ylim = c(0, max(airline, airline.bmds.d)))
#' abline(0, 1, lty = 2, col = "gray")
#' points(airline, airline.mds.d, pch = 19, col = "cyan", cex = .5)
#' points(airline, airline.bmds.d, pch = 19, col = "magenta", cex = .5)
#' legend(x = "bottomright", legend = c("Classical MDS", "Bayesian MDS"),
#'   pch = c(19, 19), col = c("cyan", "magenta"))
#' par(opar)
#'
#' # Careers of Lloyds Bank Employees
#' lloyds <- read.csv(file = system.file("extdata", "lloyds.csv",
#'   package = "dmbc"))
#' lloyds.nm <- lloyds[, 1]
#' lloyds <- lloyds[, 2:81]
#' colnames(lloyds) <- lloyds.nm
#' lloyds <- as.dist(lloyds)
#' 
#' min_p <- 1
#' max_p <- 12
#' burnin <- 200
#' nsim <- 1000
#' totiter <- burnin + nsim
#' 
#' lloyds.mds <- cmdscale(lloyds, max_p)
#' lloyds.bmds <- bmds(lloyds, min_p, max_p, burnin, nsim)
#' 
#' opar <- par(mfrow = c(1, 2))
#' plot((min_p:max_p), lloyds.bmds$mdsIC$mdsic, type = "b",
#'   main = "MDS Information Criterion", xlab = "p", ylab = "MDSIC")
#' MDSICmin <- which.min(lloyds.bmds$mdsIC$mdsic)
#' points((min_p:max_p)[MDSICmin], lloyds.bmds$mdsIC$mdsic[MDSICmin],
#'   col = "red", pch = 10, cex = 1.75, lwd = 1.5)
#' 
#' lloyds.bmds.x.mode <- bmds_get_x_mode(lloyds, lloyds.bmds, MDSICmin,
#'   min_p, max_p, start = (burnin + 1), end = totiter)
#' lloyds.bmds.d <- dist(lloyds.bmds.x.mode)
#' lloyds.mds.d <- dist(lloyds.mds[, 1:((min_p:max_p)[MDSICmin])])
#' plot(lloyds, lloyds.bmds.d, type = "n", xlab = "observed",
#'   ylab = "estimated", main = "Careers of Lloyds \n Bank Employees, 1905-1950",
#'   xlim = c(0, max(lloyds, lloyds.bmds.d)),
#'   ylim = c(0, max(lloyds, lloyds.bmds.d)))
#' abline(0, 1, lty = 2, col = "gray")
#' points(lloyds, lloyds.mds.d, pch = 19, col = "cyan", cex = .5)
#' points(lloyds, lloyds.bmds.d, pch = 19, col = "magenta", cex = .5)
#' legend(x = "topleft", legend = c("Classical MDS", "Bayesian MDS"),
#'   pch = c(19, 19), col = c("cyan", "magenta"))
#' par(opar)
#'
#' # Road distances (in km) between 21 cities in Europe
#' data(eurodist, package = "datasets")
#' 
#' min_p <- 1
#' max_p <- 10
#' burnin <- 200
#' nsim <- 1000
#' totiter <- burnin + nsim
#' 
#' eurodist.mds <- cmdscale(eurodist, max_p)
#' eurodist.bmds <- bmds(eurodist, min_p, max_p, burnin, nsim)
#' 
#' opar <- par(mfrow = c(1, 2))
#' plot((min_p:max_p), eurodist.bmds$mdsIC$mdsic, type = "b",
#'   main = "MDS Information Criterion", xlab = "p", ylab = "MDSIC")
#' MDSICmin <- which.min(eurodist.bmds$mdsIC$mdsic)
#' points((min_p:max_p)[MDSICmin], eurodist.bmds$mdsIC$mdsic[MDSICmin],
#'   col = "red", pch = 10, cex = 1.75, lwd = 1.5)
#' 
#' eurodist.bmds.x.mode <- bmds_get_x_mode(eurodist, eurodist.bmds,
#'   MDSICmin, min_p, max_p, start = (burnin + 1), end = totiter)
#' eurodist.bmds.d <- dist(eurodist.bmds.x.mode)
#' eurodist.mds.d <- dist(eurodist.mds[, 1:((min_p:max_p)[MDSICmin])])
#' plot(eurodist, eurodist.bmds.d, type = "n", xlab = "observed",
#'   ylab = "estimated", main = "Road distances (in km) \n between 21 cities in Europe",
#'   xlim = c(0, max(eurodist, eurodist.bmds.d)),
#'   ylim = c(0, max(eurodist, eurodist.bmds.d)))
#' abline(0, 1, lty = 2, col = "gray")
#' points(eurodist, eurodist.mds.d, pch = 19, col = "cyan", cex = .5)
#' points(eurodist, eurodist.bmds.d, pch = 19, col = "magenta", cex = .5)
#' legend(x = "topleft", legend = c("Classical MDS", "Bayesian MDS"),
#'   pch = c(19, 19), col = c("cyan", "magenta"))
#' par(opar)
#' }
#' @export
bmds <- function(D, min_p = 1, max_pm1 = 6, burnin = 0, nsim = 13000, ic = TRUE, verbose = TRUE) {
	if (any(is.na(D)))
		stop("NA values not allowed in the dissimilarity matrix D.")
	Dm <- as.matrix(D)
	if (is.null(n <- attr(D, "Size"))) {
		if (nrow(Dm) != ncol(Dm))
			stop("dissimilarities must be the result of the 'dist' function or a square matrix.")
	}
	
	m <- n*(n - 1)/2
	max_p <- ifelse(ic, max_pm1 + 1, max_pm1)   # an additional dimension is included when IC indexes are requested
	
	totiter <- burnin + nsim
	scale <- 2.38^2

	s_dsq <- as.numeric(sum(D^2))
	
	stress <- rmin_ssr <- numeric(max_p - min_p + 1)
	x_star <- array(NA, c((max_p - min_p + 1), n, max_p))
	x.chain <- array(NA, c(totiter, (max_p - min_p + 1), n, max_p))
	sigma2.chain <- matrix(NA, nrow = totiter, ncol = (max_p - min_p + 1))
	lambda.chain <- array(NA, c(totiter, (max_p - min_p + 1), max_p))
	accept <- matrix(NA, nrow = 2, ncol = (max_p - min_p + 1))
	for (p in min_p:max_p) {
		p.i <- p - min_p + 1
		
		if (verbose) {
			if (ic) {
				print(paste("NUMBER OF MDS DIMENSIONS: ", p, " (", p.i, "/", (max_p - min_p + 1),
          ", including an additional dimension)", sep = ""), quote = FALSE)
			} else {
				print(paste("NUMBER OF MDS DIMENSIONS: ", p, " (", p.i, "/", (max_p - min_p + 1), ")", sep = ""),
          quote = FALSE)
			}
		}
		
		accept_x_i <- accept_sigma2 <- 0
		rI_p <- diag(1, p)
		
		# initialize x
		x <- stats::cmdscale(d = D, k = p)

		# adjust the center and orientation of x
		x_adj <- adjust_x(x)
		x <- x_adj$x
		Sig_x <- x_adj$Sig_x
	
		x_int <- x
		
		# initialize sigma2 and stress
		s_res <- comp_ssr(x_int, D)
		sigma2 <- s_res/m
		stress[p.i] <- sqrt(s_res/s_dsq)
		
		# initialize lambda
		lambda <- 1/diag(Sig_x)
	
		# set parameters of the prior (use initial as prior)
		pralpha_sig <- 5
		prbeta_sig <- (pralpha_sig - 1)*sigma2
		
		alpha_lam <- 0.5
		beta_lam <- diag(Sig_x)/2
		
		pst_beta_lam <- numeric(p)
	
		# start iteration
		niter <- 0
		s_sigma2 <- 0
		sq_sigma2 <- 0
		rmin_ssr[p.i] <- stress[p.i]^2*s_dsq
		x_star[p.i, , 1:p] <- x
		
		while (niter < totiter) {
			niter <- niter + 1
			
			# generate x by using random walk Metropolis-Hastings
			for (i in 1:n) {
				cd_var <- scale*sigma2/(n - 1)*rI_p
				cd_var_inv <- 1/scale/sigma2*(n - 1)*rI_p
				cd_sig <- sqrt(scale*sigma2/(n - 1))*rI_p
				
				cdm <- x_old <- x[i, ]
				
				z <- rnorm(p)
				x_new <- z*diag(cd_sig) + cdm
				
				quad <- sum(x_old^2*lambda)
				quad_st <- sum(x_new^2*lambda)
				
				t1 <- t1_st <- s_temp <- 0
			
				for (j in 1:n) {
					if (j != i) {
						del_ij <- sqrt(sum((x_old - x[j, ])^2))
						del_st_ij <- sqrt(sum((x_new - x[j, ])^2))
						
						temp <- pnorm(del_st_ij/sqrt(sigma2), log.p = TRUE) - pnorm(del_ij/sqrt(sigma2), log.p = TRUE)
						t1 <- t1 - 0.5/sigma2*(del_ij - Dm[i, j])^2
						t1_st <- t1_st - 0.5/sigma2*(del_st_ij - Dm[i, j])^2
						s_temp <- s_temp - temp
					}
				}
				
				rl_gw <- 0
				rl_fw <- t1_st - t1 - s_temp - (quad_st - quad)/2
				
				ran_unif <- runif(1)
				if (ran_unif < exp(rl_fw - rl_gw)) {
					x_old <- x_new
					accept_x_i <- accept_x_i + 1
				}
						
				x[i, ] <- x_old
			}
		
			# adjust the center and rotation of x
			x_adj <- adjust_x(x)
			x <- x_adj$x
			Sig_x <- x_adj$Sig_x
			
			x.chain[niter, p.i, , 1:p] <- x
			
			# generate sigma2 by using random walk Metropolis-Hastings
			s_res <- comp_ssr(x, D)
			sig2_old <- sigma2
			cdvar_sig <- 2*scale*s_res^2/(((m - 2)^2)*(m - 4))
			
			ran_nor <- rnorm(1)
			sig2_new <- ran_nor*sqrt(cdvar_sig) + sig2_old
			while (sig2_new < 0) {
				ran_nor <- rnorm(1)
				sig2_new <- ran_nor*sqrt(cdvar_sig) + sig2_old
			}
			
			delta <- as.matrix(dist(x))
			rl_f <- sum(-pnorm(delta[lower.tri(delta)]/sqrt(sig2_new), log.p = TRUE) +
        pnorm(delta[lower.tri(delta)]/sqrt(sig2_old), log.p = TRUE))
			rl_f <- rl_f - (s_res/2 + prbeta_sig)*(1/sig2_new - 1/sig2_old) - (m/2 + pralpha_sig + 1)*log(sig2_new/sig2_old)
			rl_g <- pnorm(sig2_new/sqrt(cdvar_sig), log.p = TRUE) - pnorm(sig2_old/sqrt(cdvar_sig), log.p = TRUE)
			rl_w <- rl_f - rl_g
			
			ran_unif <- runif(1)
			if (ran_unif <= exp(rl_w)) {
				sig2_old <- sig2_new
				accept_sigma2 <- accept_sigma2 + 1
			}
			
			sigma2 <- sig2_old
			
			sigma2.chain[niter, p.i] <- sigma2
			
			# generate lambda using its full conditional posterior distribution
			for (ip in 1:p) {
				s1 <- sum(x[, ip]^2)
				
				pst_beta_lam[ip] <- s1/2 + beta_lam[ip]
				pst_alpha_lam <- n/2 + alpha_lam
				
				ran_gam <- rgamma(1, shape = pst_alpha_lam)
				lambda[ip] <- ran_gam/pst_beta_lam[ip]	
			}
			
			lambda.chain[niter, p.i, 1:p] <- lambda
			
			if (niter > burnin) {
				s_sigma2 <- s_sigma2 + sigma2
				sq_sigma2 <- sq_sigma2 + sigma2^2
	
				ssr <- comp_ssr(x, D)
				
				if (ssr < rmin_ssr[p.i]) {
					rmin_ssr[p.i] <- ssr
					ind_ssr <- niter
					x_star[p.i, , 1:p] <- x
				}
			}
			
			if (verbose) {
				if ((niter/100 == round(niter/100))) {
					print(paste("   iteration ", niter, "/", totiter, " ==> ",
            "acceptance x_i: ", round(accept_x_i/(n*niter), 4), " - ",
            "acceptance sigma2: ", round(accept_sigma2/niter, 4), sep = ""), quote = FALSE)
				}
			}
			
			stress[p.i] <- sqrt(rmin_ssr[p.i]/s_dsq)
			e_sigma2 <- s_sigma2/nsim
			var_sigma2 <- sq_sigma2/nsim - e_sigma2^2
		}
		accept[, p.i] <- c(accept_x_i/(n*totiter), accept_sigma2/totiter)
	}
	
	# compute MDSIC
	if (ic) {
		mdsIC <- mdsic(x_star, rmin_ssr, n, min_p, max_p)
	} else {
		mdsIC <- NA
	}
	
	return(list(x.chain = x.chain, sigma2.chain = sigma2.chain, lambda.chain = lambda.chain,
              stress = stress, mdsIC = mdsIC, accept = accept))
}

#' Information criterion for Bayesian multidimensional scaling (BMDS).
#' 
#' \code{mdsic} computes the information criterion for a set of Bayesian
#'   multidimensional scaling (BMDS) solutions using the approach in
#'   Oh & Raftery (2001).
#' 
#' @param x_star An array containing the latent configurations
#'   estimated using \code{\link{bmds}}.
#' @param rmin_ssr A numeric vector providing the ratios of SSR
#'   for the latent dimensions requested.
#' @param n A length-one numeric vector providing the number of objects.
#' @param min_p A length-one numeric vector providing the minimum value
#'   of the latent space dimension to use.
#' @param max_p A length-one numeric vector providing the maximum
#'   value of the latent space dimension to use.
#' @return A list with the following elements:
#'   \describe{
#'     \item{\code{mdsic}}{A numeric vector with the values of MDSIC index.}
#'     \item{\code{bic}}{A numeric vector with the values of the BIC index.}
#'   }
#' @author Sergio Venturini \email{sergio.venturini@unito.it}
#' @seealso \code{\link{bmds}} for Bayesian (metric) multidimensional scaling
#'   and \code{\link{comp_ssr}} for the computation of SSR.
#' @references
#'   Oh, M.-S., Raftery, A. E. (2001), "Bayesian Multidimensional Scaling and
#'   Choice of Dimension", Journal of the American Statistical Association,
#'   96, 1031-1044.
#' @examples
#' \dontrun{
#' # Road distances (in km) between 21 cities in Europe
#' data(eurodist, package = "datasets")
#' 
#' min_p <- 1
#' max_p <- 10
#' burnin <- 200
#' nsim <- 1000
#' totiter <- burnin + nsim
#' 
#' eurodist.mds <- cmdscale(eurodist, max_p)
#' eurodist.bmds <- bmds(eurodist, min_p, max_p, burnin, nsim)
#' 
#' plot((min_p:max_p), eurodist.bmds$mdsIC$mdsic, type = "b",
#'   main = "MDS Information Criterion", xlab = "p", ylab = "MDSIC")
#' MDSICmin <- which.min(eurodist.bmds$mdsIC$mdsic)
#' points((min_p:max_p)[MDSICmin], eurodist.bmds$mdsIC$mdsic[MDSICmin],
#'   col = "red", pch = 10, cex = 1.75, lwd = 1.5)
#' }
#' @export
mdsic <- function(x_star, rmin_ssr, n, min_p = 1, max_p = 6) {
	Rj <- MDSIC <- BIC <- BIC_sum <- numeric(max_p - min_p)
	m <- n*(n - 1)/2

	for (p in min_p:(max_p -1)) {
		p.i <- p - min_p + 1

		if (p.i == 1) rLL1 <- -m*log(rmin_ssr[p.i])/2
		
		rLRT <- (m - 2)*(log(rmin_ssr[p.i + 1]) - log(rmin_ssr[p.i]))
		
		s_j_pp1 <- apply(drop(x_star[p.i + 1, , ]^2), 2, sum)
		s_j_p <- apply(drop(x_star[p.i, , ]^2), 2, sum)
		rs_j <- s_j_pp1/s_j_p
		s <- sum(log(rs_j*(n + 1)/(n + rs_j)), na.rm = TRUE)
		penalty <- (n + 1)*(s + log(n + 1))	 
		
		Rj[p.i] <- rLRT + penalty
		BIC[p.i] <- rLRT - (n + 1)*log(m)/2
		if (p.i == 1) {
			MDSIC[p.i] <- (m - 2)*log(rmin_ssr[1])
			BIC_sum[p.i] <- -2*rLL1 + (n + 1)*log(m)
		} else {
			MDSIC[p.i] <- MDSIC[1] + sum(Rj[1:(p.i - 1)])
			BIC_sum[p.i] <- -2*(sum(BIC[1:(p.i - 1)]) + rLL1)
		}
	}
	
	return(list(mdsic = MDSIC, bic = BIC_sum))
}
