#' @title Maximum likelihood rate estimation for traits and phylogenies 
#' @description Function for the maximum likelihood estimation of rate parameters on a trait and phylogeny.
#' @param rateData an object of class \code{rateData}
#' @param rate a vector of relative rate parameters. The length of the vector is equal to the number of rates being estimated. If \code{rate=NULL} then rates are equal.
#' @param fixed A vector stating whether each parameter should be allowed to vary (either \code{FALSE} which results in a start value of 1, or a numeric start value) or should be fixed (\code{TRUE}).
#' @param rateMIN Minimum value for the rate parameters
#' @param rateMAX Maximum value for the rate parameters
#' @param common.mean a logical specififying whether each rate category should have its own mean (\code{common.mean=FALSE}) or all categories should have the same mean (\code{common.mean=FALSE}). See Thomas et al. (2009) for a discussion on the impact of assumptions about mean on rate estimates.
#' @param lambda.est Logical. Fit Pagel's lambda.
#' @param meserr Logical. Include measurement error.
#' @return MLRate Maximum likelihood estimates of the rate parameters
#' @return Max.lik  Maximum (log) likeihood
#' @return AIC AIC for maximum likelihood model
#' @return AICc  AICc for maximum likelihood model
#' @return convergence convergence value from \code{optim}
#' @return n.parameters Number of parameters in the model (how many means and rate categories)
#' @references Thomas GH, Freckleton RP, & Szekely T. 2006. Comparative analyses of the influence of developmental mode on phenotypic diversification rates in shorebirds. Proceedings of the Royal Society B 273, 1619-1624.
#' @references Thomas GH, Meiri S, & Phillimore AB. 2009. Body size diversification in Anolis: novel environments and island effects. Evolution 63, 2017-2030.
#' @author Gavin Thomas
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
#' # A model with a different rate in each of the four groups. The 'fixed' command is used to determine
#' # whether a particular rate is to be constrained or not. Use '1' to fix a group and 'FALSE' to show
#' # that the parameter is not fixed and should be estimated. The values should be entered in the same 
#' # order as the ranking of the groups. That is, group 0 (small islands) takes position one in the 
#' # fixed vector, group 1 (large island trunk crown and trunk ground) takes position 2 and so on. 
#' # The default is to allow each group to take a different mean. 
#'
#' optim.likRatePhylo(anolis.rateData, rate=c(1,1,1,1), common.mean=TRUE, lambda.est=FALSE)
#' @export

optim.likRatePhylo <-
function(rateData, rate=NULL, fixed = NULL, rateMIN = 0.001, rateMAX = 50, common.mean=FALSE, lambda.est=TRUE, meserr=FALSE) {
		
		if(is.null(rate))  { rate <- c(rep(1,length(rateData$Vmat))) } else { rate <- rate }	
	
		if(is.null(fixed))  { op <- fixed <- c(TRUE, rep(FALSE,length(rateData$Vmat) - 1)) } else { op <- fixed }
				
		lower <- c(rep(rateMIN, length(rateData$Vmat)), 0.00001)
		upper <- c(rep(rateMAX, length(rateData$Vmat)), 1)
	
		mvl <- make.likRatePhylo(rateData, fixed=op, common.mean=common.mean, lambda.est=lambda.est, meserr=meserr)

			vo <- try(optim(c(rate, 1, 1), mvl, method = "L-BFGS-B", lower = lower, upper = upper))
			MLRate <- vo$par[1:length(fixed)]
			Lambda <- vo$par[1+length(fixed)]
			
			fixed[which(fixed==FALSE)] <- MLRate[which(fixed==FALSE)]
			MLRate <- fixed

			ML <- -vo$value
			convergence <- vo$convergence
			n <- length(rateData$y)
			
				if(length(op)!=length(which(op==FALSE))) {
					if(common.mean==TRUE) {k <- 2 + length(which(op==FALSE) + lambda.est + meserr)
							} else { k <- (length(which(op==FALSE)) +1 + length(op)) + lambda.est + meserr }
							
							} else {
					
					if(common.mean==TRUE) {k <- 1 + length(which(op==FALSE) + lambda.est + meserr)
							} else { k <- (length(which(op==FALSE)) + length(op)) + lambda.est + meserr}
							}
							
			aic <- -2 * ML + 2 * k
			aicc <- -2 * ML + 2 * k + ((2*k*(k+1))/(n-k-1))
			
			
	ML.RatePhylo <- list(MLRate = MLRate, Lambda = Lambda, Max.lik = ML, aic = aic, aicc = aicc, convergence=convergence, n.parameters = k)
	return(ML.RatePhylo)
	
}

