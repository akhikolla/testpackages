#' @export glo
glo <- texmexFamily(name = 'GLO',
					param = c(mu=0, phi=0, xi=0),
					log.lik = function(data, ...) {
						y <- data$y
						X.mu <- data$D$mu
						X.phi <- data$D$phi
						X.xi <- data$D$xi

						n.mu <- ncol(X.mu)
						n.phi <- n.mu + ncol(X.phi)
						n.end <- n.phi + ncol(X.xi)

						function(param) {
							stopifnot(length(param) == n.end)
							mu <- X.mu %*% param[1:n.mu]
							phi <- X.phi %*% param[(1 + n.mu):n.phi]
							xi <- X.xi %*% param[(1 + n.phi):n.end]
							if (any(1 + xi/exp(phi)*(y-mu) <= 0)) -Inf
							else sum(dglo(y, mu, exp(phi), xi, log.d=TRUE))
						}
					}, # Close log.lik
					info = NULL, # will mean that numerical approx gets used
					sandwich = NULL,
					delta = function(param, m, model){ # model not used but required by a calling function
						out <- rep(1, 3)
						out[2] <- exp(param[2])/param[3] * ((m-1)^(-param[3])  - 1) # Coles p.56
						out[3] <- -exp(param[2]) * param[3]^(-2) * ((m-1)^(-param[3])  - 1) + exp(param[2])/param[3] * (m-1)^(-param[3]) * log(m-1)
						out
					}, # Close delta
					start = function(data){
						y <- data$y
						X.mu <- data$D[[1]]
						X.phi <- data$D[[2]]
						X.xi <- data$D[[3]]

						c(mean(y), rep(0, ncol(X.mu)-1), log(IQR(y)/2),
						  rep(.001, -1 + ncol(X.phi) + ncol(X.xi)))
					}, # Close start
					endpoint = function(param, model){
						ep <- param[, 1] - exp(param[, 2]) / param[, 3]
						ep[!(param[,3]<0)] <- Inf
						ep
					},
					rng = function(n, param, model){
						rglo(n, c(param[, 1]), exp(c(param[, 2])), c(param[, 3]))
					},
					density = function(x, param, model, log.d=FALSE){
						dglo(x, c(param[, 1]), exp(c(param[, 2])), c(param[, 3]), log.d=log.d)
					},
					prob = function(x, param, model){
						pglo(x, c(param[, 1]), exp(c(param[, 2])), c(param[, 3]))
					},
					quant = function(p, param, model){
						qglo(p, c(param[, 1]), exp(c(param[, 2])), c(param[, 3]))
					},
					resid = function(o) { # these have a standard Logistic distribution under the model
						p <- texmexMakeParams(coef(o), o$data$D)
						shift <- (o$data$y - p[,1]) / exp(p[,2])
						standard.logistic <- .log1prel(p[,3]*shift)*shift
						standard.logistic
					}, # Close resid
					rl = function(m, param, model){
						qglo(1/m, param[,1], exp(param[,2]), param[,3], lower.tail=FALSE)
					}
)
#' Generalized logistic distribution
#'
#' @description Density, distribution and quantile functions, and random number
#'   generation for the Generalized logistic distribution
#' @param x,q,p Value, quantile or probability respectively.
#' @param n Number of random numbers to generate.
#' @param mu Location parameter.
#' @param sigma Scale parameter.
#' @param xi Shape parameter.
#' @param log.d,log.p Whether to work on the log scale.
#' @param lower.tail Whether to return the lower tail.
#' @family pglo qglo dglo rglo
#' @aliases pglo qglo dglo rglo
#' @export
#' @rdname dglo
rglo <- function(n, mu, sigma, xi){
	## use standard GLO ~ exp(xi*S - 1)/xi
	## where S is a standard Logistic

	s <- rlogis(n)

	## expand mu, sigma, and xi to be n long
	## this is necessary to ensure that we get
	## exactly n random numbers if mu, sigma, xi
	## are greater than n long
	n     <- length(s)
	mu    <- rep(mu, length.out=n)
	sigma <- rep(sigma, length.out=n)
	xi    <- rep(xi, length.out=n)

	## and here we go
	standard.glo <- .exprel(xi*s)*s

	if( sum(xi == 0)){
		standard.glo[xi==0] <- s[xi==0]
	}

	mu + sigma * standard.glo
}

#' @export
#' @name dglo
dglo <- function(x, mu, sigma, xi, log.d=FALSE){
	## shift and scale
	x <- (x - mu) / sigma

	logrel <- .log1prel(xi*x) * x #Accurately compute log(1 + x) / x

	log.density <- -log(sigma) - log1p(xi*x) - logrel - 2* log1p(exp(-logrel))

	## make exp(Inf) > Inf
	log.density[logrel==(-Inf)] <- -Inf

	if (!log.d) {
		exp(log.density)
	} else {
		log.density
	}
}

#' @export
#' @rdname dglo
pglo <- function(q, mu, sigma, xi, lower.tail=TRUE, log.p=FALSE){
		## first shift and scale
		q <- (q - mu) / sigma

		## now set the lengths right
		n  <- max(length(q), length(xi))
		q  <- rep(q, length.out=n)
		xi <- rep(xi, length.out=n)

		res <- .log1prel(xi*q) * q

		logP <- -log1p(exp(res))

		if(log.p){
			if(lower.tail) return(log1p(-exp(logP)))else return(logP)
		} else {
			if(lower.tail) return(1-exp(logP)) else return(exp(logP))
		}
}

#' @export
#' @rdname dglo
qglo <- function(p, mu, sigma, xi, lower.tail=TRUE, log.p=FALSE){

	if ((!log.p) && any((p < 0) || (p > 1))) {
		stop("p must lie between 0 and 1 if log.p=FALSE")
	}

	if(!lower.tail) {
		if(!log.p) {
			p <- 1-p
		} else {
			p <- log(1-exp(p))
		}
	}

	if(log.p) {
		p <- exp(p)
	}

	standard <- qlogis(p=p)

	if( sum(xi != 0)){
		standard[xi != 0] <- 1/xi * (1/(1/p-1)^xi - 1)
	}

	## and now shift and scale
	mu + sigma * standard
}


