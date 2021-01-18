#' @title Maximum likelihood rate estimation for traits and phylogenies
#' @description Full function for maximum likelihood estimation of rate parameters and comparison to a single rate model.
#' @param rateData an object of class \code{rateData}
#' @param rate a vector of relative rate parameters. The length of the vector is equal to the number of rates being estimated. If \code{rate=NULL} then rates are equal.
#' @param fixed A vector stating whether each parameter should be allowed to vary (either \code{FALSE} which results in a start value of 1, or a numeric start value) or should be fixed (\code{TRUE}).
#' @param pretty Display the output nicely (\code{pretty=TRUE}) or as a list (\code{pretty=FALSE})
#' @param rateMIN Minimum value for the rate parameters.
#' @param rateMAX Maximum value for the rate parameters
#' @param common.mean a logical specififying whether each rate category should have its own mean (\code{common.mean=FALSE}) or all categories should have the same mean (\code{common.mean=TRUE}). See Thomas et al. (2009) for a discussion on the impact of assumptions about mean on rate estimates.
#' @param lambda.est Logical. Estimate Pagel's lambda.
#' @param est.CI Logical. Estimate approximate confidence intervals for rate parameters.
#' @param meserr Logical. Incorporate measurement error.
#' @param file File string for output. Only used if \code{pretty=TRUE}.
#' @return If \code{pretty=FALSE}, returns a list containing:
#' \itemize{
#' \item {MLRate} Maximum likelihood estimates of the rate parameters
#' \item {Lambda}  Maximum likelihood estimate of lambda
#' \item {LCI} Approximate lower confidence intervals for rate
#' \item {UCI} Approximate upper confidence intervals for rate parameters
#' \item {means} Means for each category
#' \item {nParam} Number of parameters in the model (how many means and rate categories)
#' \item {Max.lik} Maximum (log) likeihood
#' \item {AIC} for maximum likelihood model
#' \item {AICc} for maximum likelihood model
#' \item {LambdaSingle} Maximum likelihood estimate of lambda for the single rate model
#' \item {Lik1}  Likelihood of the equivalent single rate model
#' \item {Likelihood} ratio statistic of "Max.lik" vs "Lik1"
#' \item {P}  P values for the LR statistic
#' \item {df} Degrees of freedom for the LR statistic
#' \item {AIC.rate1} AIC for single rate model
#' \item {AICc.rate1} AICc for single rate model
#' }
#' @return If \code{pretty=TRUE}, prints a nice version of the list to screen. If \code{file} is specified the pretty output will be sent to file, not the console.
#' @note Unlike phyloMean and likRatePhylo (that use treatment contrasts), the means reported here are the actual values
#' @references Thomas GH, Freckleton RP, & Szekely T. 2006. Comparative analyses of the influence of developmental mode on phenotypic diversification rates in shorebirds. Proceedings of the Royal Society B 273, 1619-1624.
#' @references Thomas GH, Meiri S, & Phillimore AB. 2009. Body size diversification in Anolis: novel environments and island effects. Evolution 63, 2017-2030.
#' @author Gavin Thomas, Rob Freckleton
#' @examples ## Read in phylogeny and data from Thomas et al. (2009)
#' @examples data(anolis.tree)
#' @examples data(anolis.data)
#' @examples ## Convert data to class rateData with a rateMatrix object as input
#' @examples anolis.rateMatrix <- as.rateMatrix(phy=anolis.tree, x="geo_ecomorph",
#' @examples data=anolis.data)
#' @examples anolis.rateData <- as.rateData(y="Female_SVL", x="geo_ecomorph", 
#' @examples rateMatrix = anolis.rateMatrix, phy=NULL, data=anolis.data, log.y=TRUE)
#' @examples # A model with a different rate in one of the four groups. The 'fixed' command is used to determine
#' @examples # whether a particular rate is to be constrained or not. Use '1' to fix a group and 'FALSE' to show
#' @examples # that the parameter is not fixed and should be estimated. The values should be entered in the same 
#' @examples # order as the ranking of the groups. That is, group 0 (small islands) takes position one in the 
#' @examples # fixed vector, group 1 (large island trunk crown and trunk ground) takes position 2 and so on. 
#' @examples # The default is to allow each group to take a different mean. 
#' @examples ML.RatePhylo(anolis.rateData, fixed=c(1, FALSE, FALSE, FALSE), pretty=TRUE, lambda.est = FALSE)
#' @export

ML.RatePhylo <-
function(rateData, rate=NULL, fixed = NULL, pretty = TRUE, rateMIN = 0.001, rateMAX = 50, common.mean=FALSE, lambda.est=TRUE, est.CI=FALSE, meserr=FALSE, file=NULL) {
						
		if(is.null(rate))  { rate <- c(rep(1,length(rateData$Vmat))) } else { rate <- rate }	
	
		if(is.null(fixed))  { op <- fixed <- c(rep(FALSE,length(rateData$Vmat) - 1), TRUE) } else { op <- fixed }		
		
				ovp <- optim.likRatePhylo(rateData, rate, fixed, rateMIN, rateMAX, common.mean=common.mean, lambda.est=lambda.est, meserr=meserr)
				max.lik.rate <- ovp$MLRate
				lambda <- ovp$Lambda
				max.lik <- ovp$Max.lik
				singlerate <- optim.likRatePhylo(rateData, rate, rep(TRUE, length(rateData$Vmat)), rateMIN, rateMAX, common.mean=common.mean, lambda.est=lambda.est, meserr=meserr)
				lik1 <- singlerate$Max.lik
				lambda_single <- singlerate$Lambda
				D <- 2 * (max.lik - lik1)
	
	
				mlparams <- likRatePhylo(rateData, rate=max.lik.rate, common.mean=common.mean, lambda.est=lambda)
	
	
				if (common.mean==TRUE) { mu <- mlparams$mu } else {
					
					mu <- rep(NA, length(mlparams$mu))
					mu[1] <- mlparams$mu[1]
					for (i in 2:length(mlparams$mu)) { mu[i] <- mlparams$mu[1] + mlparams$mu[i] }
				}
		
														
		
				
				if (length(op) == length(which(op==FALSE)))  { k2 <- length(op) - 1 
						} else { k2 <- length(which(op==FALSE)) }
				
				if(length(op)!=length(which(op==FALSE))) {
					if(common.mean==TRUE) {k <- 2 + length(which(op==FALSE)) + lambda.est + meserr
							} else { k <- (length(which(op==FALSE)) +1 + length(op)) + lambda.est + meserr}
							
							} else {
					
					if(common.mean==TRUE) {k <- 1 + length(which(op==FALSE)) + lambda.est + meserr
							} else { k <- (length(which(op==FALSE)) + length(op)) + lambda.est + meserr}
							}
								
				pval <- 1- pchisq(D, k2)
				n <- length(rateData$y)
				
				
	if (est.CI) { CIs.rate <- RatePhylo.allCI(rateData, max.lik.rate, fixed=rep("FALSE", length(rateData$Vmat)), common.mean=common.mean)
	} else { CIs.rate <- matrix(NA, ncol=2, nrow=1) }
				aic <- -2 * max.lik + 2 * k
				aicc <- -2 * max.lik + 2 * k + ((2*k*(k+1))/(n-k-1))
				
				if(common.mean==TRUE) {
						aic.rate1 <- -2 * lik1 + 2 * 2
						aicc.rate1 <- -2 * lik1 + 2 * 2 + ((2*2*(2+1))/(n-2-1))
					} else {
						aic.rate1 <- -2 * lik1 + 2 * (1 + length(op))
						aicc.rate1 <- -2 * lik1 + 2 * (1 + length(op)) + ((2*(1 + length(op))*((1 + length(op))+1))/(n-(1 + length(op))-1))
						}
				
				if (is.null(file)) { file <- "" }
				if(pretty == TRUE) {
					cat("____________________________\n", file=paste(file), append=TRUE)
					cat("Maximum likelihood estimation: rates:\n\n", file=paste(file), append=TRUE)
					cat("Lambda:          ", lambda, "\n", file=paste(file), append=TRUE)
					cat("Brownian variance (rate): ", mlparams$s2, "\n", file=paste(file), append=TRUE)
					cat("ML estimates of group relative rates :          ", max.lik.rate, "\n", file=paste(file), append=TRUE)
					cat("Lower confidence intervals for rates:", CIs.rate[,1], "\n", file=paste(file), append=TRUE)
					cat("Upper confidence intervals for rates:", CIs.rate[,2], "\n", file=paste(file), append=TRUE)
					cat("ML estimates of group means : ", mu, "\n\n", file=paste(file), append=TRUE)
					cat("Number of parameters: ", k, "\n", file=paste(file), append=TRUE)
					cat("Maximised log likelihood: ", max.lik, "\n", file=paste(file), append=TRUE)
					cat("  AIC = ", aic, "  \n", file=paste(file), append=TRUE)
					cat("  AICc = ", aicc, "  \n\n", file=paste(file), append=TRUE)
					cat("____________________________\n", file=paste(file), append=TRUE)
					cat("Comparison with single rate model\n", file=paste(file), append=TRUE)
					cat("Lambda (single rate model):          ", lambda_single, "\n", file=paste(file), append=TRUE)
					cat("Log likelihood (single rate): ", lik1, "\n", file=paste(file), append=TRUE)
					cat("LR statistic (ML rates vs single rate):", D, file=paste(file), append=TRUE)
					cat("  P = ", pval , file=paste(file), append=TRUE)
					cat(" df = ", k2, " \n", file=paste(file), append=TRUE)
					cat("  Single rate AIC = ", aic.rate1, "  \n", file=paste(file), append=TRUE)
					cat("  Single rate AICc = ", aicc.rate1, "  \n", file=paste(file), append=TRUE)
					cat("____________________________\n", file=paste(file), append=TRUE)
					}
				if(pretty == FALSE) {
					max.lik.rate.list <- list(BRvar=mlparams$s2, MLRate = max.lik.rate, LCI = CIs.rate[,1], UCI = CIs.rate[,2],  means=mu, nParam = k, lambda=lambda, Max.lik = max.lik, AIC = aic, AICc=aicc, LambdaSingle = lambda_single, Lik1 = lik1, LR = D, P = pval, df = k2, AIC.rate1=aic.rate1, AICc.rate1=aicc.rate1)
					return(max.lik.rate.list) }	
	}

