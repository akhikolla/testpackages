geom_phi_log_prob <- function(phi, p, using_terms) {
	x <- length(phi)
	(x-1) * log(1-p) + log(p) - lchoose(length(using_terms), x)
}

#' @importFrom stats dgamma
discrete_gamma <- function(using_terms, max_size=10, shape=10, rate=4) {
	n_terms <- length(using_terms)
	dist <- dgamma(seq(max_size), shape=shape, rate=rate)
	dens <- log(dist/sum(dist))

	function(phi) {
		x <- length(phi)
		if (x > max_size) -Inf
		else dens[x] - lchoose(n_terms, x)
	}
}


