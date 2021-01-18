#' @title Conversion among data and phylogeny objects
#' @description Transforms the expected variance and covariances among species according to hypotheses of rate variation between lineages.
#' @param rateData an object of class \code{rateData}
#' @param rate a vector of relative rate parameters. The length of the vector is equal to the number of rates being estimated. If \code{rate=NULL} then rates are equal.
#' @return retMat Rate-transformed variance covariance matrix
#' @references Thomas GH, Freckleton RP, & Szekely T. 2006. Comparative analyses of the influence of developmental mode on phenotypic diversification rates in shorebirds. Proceedings of the Royal Society B 273, 1619-1624.
#' @references Thomas GH, Meiri S, & Phillimore AB. 2009. Body size diversification in Anolis: novel environments and island effects. Evolution 63, 2017-2030.
#' @author Gavin Thomas
#' @examples
#'  data(anolis.tree)
#'  data(anolis.data)
#'  ## Convert data to class rateData with a rateMatrix object as input
#'  anolis.rateMatrix <- as.rateMatrix(phy=anolis.tree, x="geo_ecomorph", data=anolis.data)
#'  anolis.rateData <- as.rateData(y="Female_SVL", x="geo_ecomorph", 
#'  rateMatrix = anolis.rateMatrix, phy=NULL, data=anolis.data, log.y=TRUE)
#'  
#'  # Tranform the expected variance covariance matrix so that the rates in the first and last 
#'  # categories are equal (both 1) whereas the rate in the second category is twice as fast (2) and 
#'  # the rate in the third category is ten times slower.
#'  
#'  trans.anolis.rateData <- transformRateMatrix(rateData=anolis.rateData, rate = c(1,2,0.1,1))
#' @export


transformRateMatrix <- function(rateData, rate=NULL) {
			V <- rateData$Vmat
			
			if(is.null(rate))  { rate <- c(rep(1,length(rateData$Vmat))) } else { rate <- rate }	
			
			nV <- length(rate)
			
			if(length(rate) != length(rateData$Vmat)){stop("The number of rates defined differs from the number of rate matrices")}
			
			v1 <- V[[1]]		
			rateMats <- vector(mode="list", length = nV)
			retMat <- matrix(0, nrow = dim(v1)[1], ncol = dim(v1)[2])
			
			for(i in 1:nV) {
			   rateMats[[i]] <- rate[i] * V[[i]]  
			   retMat <- retMat + rateMats[[i]]
					}
	
	return(retMat)}
