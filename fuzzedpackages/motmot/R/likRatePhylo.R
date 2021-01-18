#' Log-likelihood rate estimation for traits and phylogenies
#'
#' This function calculates the log-likelihood, phylogenetic mean, and Brownian variance for a trait and a phylogeny transformed according to variation in relative rates. 
#' @param rateData an object of class \code{rateData}
#' @param rate a vector of relative rate parameters. The length of the vector is equal to the number of rates being estimated. If \code{rate=NULL} then rates are equal.
#' @param common.mean a logical specififying whether each rate category should have its own mean (\code{common.mean=FALSE}) or all categories should have the same mean (\code{common.mean=FALSE}). See Thomas et al. (2009) for a discussion on the impact of assumptions about mean on rate estimates.#' 
#' @param lambda.est Logical. Fit Pagel's lambda.
#' @param lambda Logical. Numeric value for lambda from 0-1.
#' @param meserr Logical. Logical. Include measurement error.
#' @param sigmaScale Logical. Scalar for measurement error relative to tree.
#' @return ll log-likelihood of the model
#' @return mu phylogenetically corrected mean(s)
#' @return s2 Brownian variance
#' @references Thomas GH, Freckleton RP, & Szekely T. 2006. Comparative analyses of the influence of developmental mode on phenotypic diversification rates in shorebirds. Proceedings of the Royal Society B 273, 1619-1624.
#' @references Thomas GH, Meiri S, & Phillimore AB. 2009. Body size diversification in Anolis: novel environments and island effects. Evolution 63, 2017-2030.
#' @author Gavin Thomas, Rob Freckleton
#' @note The means are output as treatment contrasts.
#' @examples
#' data(anolis.tree)
#' data(anolis.data)
#' 
#' ## Convert data to class rateData with a rateMatrix object as input
#' anolis.rateMatrix <- as.rateMatrix(phy=anolis.tree, x="geo_ecomorph", data=anolis.data)
#' 
#' anolis.rateData <- as.rateData(y="Female_SVL", x="geo_ecomorph", 
#' rateMatrix = anolis.rateMatrix, phy=NULL, data=anolis.data, log.y=TRUE)
#' 
#' ## Calculate phylogenetic mean, variance, log likelihood for a model where the first 
# rate category (small islands) takes position one in the rate vector, group 1
# (large island trunk crown and trunk ground) takes position 2 and so on.
# The rates in the first and last categories are equal (both 1) whereas the rate in the
# second category is twice as fast (2) and the rate in the third category is ten times slower.
# Means are allowed to differ. 
#' 
#' # mean only
#' phyloMean(rateData=anolis.rateData, rate = c(1,2,0.1,1), common.mean = FALSE)
#' 
#' # variance only
#' phyloVar(rateData=anolis.rateData, rate = c(1,2,0.1,1), common.mean = FALSE) 
#' 
#' # mean, variance and log-likelihood
#' likRatePhylo(rateData=anolis.rateData, rate = c(1,2,0.1,1), common.mean = FALSE) 
#' @export

likRatePhylo <-
function(rateData, rate=NULL, common.mean=FALSE, lambda.est=TRUE, lambda=1, meserr=FALSE, sigmaScale=NULL) {
	
		if(is.null(rate))  { rate <- c(rep(1,length(rateData$Vmat))) } else { rate <- rate }	

			y <- rateData$y
			x <- as.factor(rateData$x)
			if (common.mean==FALSE) {k <- nlevels(x)} else {k <- 1}
	
			V <- transformRateMatrix(rateData, rate)

			x <- make.anc(y, x)
						
			if (!lambda.est & !meserr) {
				logDetV <- determinant(V)$modulus
				mu <- phyloMean(rateData, rate, common.mean, lambda.est, lambda)
				s2 <- phyloVar(rateData, rate, common.mean, lambda.est, lambda)
				n <- length(x[,1])
				ll <- -n / 2.0 * log( 2 * pi) - n / 2.0 * log(s2) - logDetV / 2.0 - (n - k)/2.0
				lik.RatePhylo <- ( list(ll = ll, mu = mu, s2 = s2) )
			}

			if (lambda.est & !meserr) {
				v.temp <- V
				diag(v.temp) <- rep(0, dim(V)[1])
				V.lam <- lambda*v.temp
				diag(V.lam) <- diag(V)
				V <- V.lam
				logDetV <- determinant(V)$modulus
				mu <- phyloMean(rateData, rate, common.mean, lambda.est, lambda)
				s2 <- phyloVar(rateData, rate, common.mean, lambda.est, lambda)
				n <- length(x[,1])
				ll <- -n / 2.0 * log( 2 * pi) - n / 2.0 * log(s2) - logDetV / 2.0 - (n - k)/2.0
				lik.RatePhylo <- ( list(ll = ll, mu = mu, s2 = s2) )
			}

	
			if (!lambda.est & meserr) {
				
				if (is.null(sigmaScale)) { stop("Estimate of sigmaScale required") }
				
				if(common.mean==FALSE) {x <- x} else { x <- rep(1, length(x[,1]))}

				diag(V) <- diag(V) + rateData$meserr/sigmaScale
				logDetV <- determinant(V)$modulus
				
				iV <- solve(V)
				xVix <- crossprod(x, iV %*% x)
				xViy <- crossprod(x, iV %*% y)
				mu <- solve(xVix) %*% xViy 

				e <- y - x %*% mu
				s2 <- crossprod(e, iV %*% e)
				n <- length(y) 
				phylo.var <- ( s2 / (n - k) )
				
				
				ll <- -n / 2.0 * log( 2 * pi) - n / 2.0 * log(phylo.var) - logDetV / 2.0 - (n - k)/2.0
				lik.RatePhylo <- ( list(ll = ll, mu = mu, s2 = phylo.var) )
			}
	
	
			if (lambda.est & meserr) {
		
				if (is.null(sigmaScale)) { stop("Estimate of sigmaScale required") }
		
				if(common.mean==FALSE) {x <- x} else { x <- rep(1, length(x[,1]))}
		
				
				v.temp <- V
				diag(v.temp) <- rep(0, dim(V)[1])
				V.lam <- lambda*v.temp
				diag(V.lam) <- diag(V)
				V <- V.lam
				diag(V) <- diag(V) + rateData$meserr/sigmaScale
				logDetV <- determinant(V)$modulus
		
				iV <- solve(V)
				xVix <- crossprod(x, iV %*% x)
				xViy <- crossprod(x, iV %*% y)
				mu <- solve(xVix) %*% xViy 
		
				e <- y - x %*% mu
				s2 <- crossprod(e, iV %*% e)
				n <- length(y) 
				phylo.var <- ( s2 / (n - k) )
		
		
				ll <- -n / 2.0 * log( 2 * pi) - n / 2.0 * log(phylo.var) - logDetV / 2.0 - (n - k)/2.0
				lik.RatePhylo <- ( list(ll = ll, mu = mu, s2 = phylo.var, lambda=lambda) )
			}
	
	
				return(lik.RatePhylo)
	}






