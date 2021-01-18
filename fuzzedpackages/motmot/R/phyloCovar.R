#' @title Calculation of Brownian (co)variance using independent contrasts.
#' @description Calculates the Brownian variance (single trait) or variance-covariance matrix (mutliple traits) using phylogenetically independent contrasts.
#' @param x A continuous trait
#' @param phy An object of class \code{phylo} (see \pkg{ape}).
#' @param estimator Should Brownian variance (or covariance) be based on the unbiased ("unbiased" - default) or maximum likelihood ("ML") estimator.
#' @return brownianVariance Brownian variance (or covariance for multiple traits) given the data and phylogeny
#' @references Felsenstein J. 1973. Maximum-likelihood estimation of evolutionary trees from continuous characters. Am. J. Hum. Genet. 25, 471-492.
#' @references Felsenstein J. 1985. Phylogenies and the comparative method. American Naturalist 125, 1-15.
#' @references Freckleton RP & Jetz W. 2009. Space versus phylogeny: disentangling phylogenetic and spatial signals in comparative data. Proc. Roy. Soc. B 276, 21-30.
#' @author Gavin Thomas, Rob Freckleton
#' @examples 
#' data(anolis.tree)
#' data(anolis.data)
#' ## calculate Brownian variance of female SVL
#' female.svl <- matrix(anolis.data[,"Female_SVL"], 
#' dimnames=list(rownames(anolis.data)))
#' input.data <- sortTraitData(phy=anolis.tree, y=female.svl, log.trait=TRUE)
#' phyloCovar(x=input.data$trait, phy = input.data$phy)
#' @export

phyloCovar <- function(x, phy, estimator="unbiased") {
		
		if (is.matrix(x)==FALSE) { stop("Trait data must be a matrix with taxon names as row names")}
		
		n <- length(phy$tip.label)
		phy <- reorder(phy, order = "pruningwise")
	
		x <- as.matrix(x[phy$tip.label,])
		
		contrasts <- apply(x, MARGIN=2, FUN=pic.motmot, phy=phy)
		rawVariances <- c(contrasts[[1]]$contr[,2], contrasts[[1]]$V)
		rawContrasts <- matrix(NA, nrow=n, ncol=ncol(x))
		
		for (i in 1:ncol(x)) {
			rawContrasts[,i] <- c(contrasts[[i]]$contr[,1],0)
			}
			
		brCov <- matrix(NA, nrow=ncol(x), ncol=ncol(x))
		
	if (estimator=="unbiased") {
		for (i in 1:ncol(x)) {
			for (k in 1:ncol(x)) {	
				brCov[k,i] <- brCov[i,k] <- crossprod(rawContrasts[,k]/sqrt(rawVariances), rawContrasts[,i]/sqrt(rawVariances)) / (n-1)
								}}
	}
	
	if (estimator=="ML") {
		for (i in 1:ncol(x)) {
			for (k in 1:ncol(x)) {	
				brCov[k,i] <- brCov[i,k] <- crossprod(rawContrasts[,k]/sqrt(rawVariances), rawContrasts[,i]/sqrt(rawVariances)) / (n)
			}}
	}
		
		return(brCov)
		}
		

