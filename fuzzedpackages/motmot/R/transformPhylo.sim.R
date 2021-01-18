#' @title Phylogenetic tree transformations
#' @description Simulates trait data on a tree using a specified model of evolution (see details).
#' @param phy An object of class \code{phylo} (see \pkg{ape}).
#' @param n Number of simulations
#' @param x Vector, matrix or data.frame (with taxon names as names or rownames) of categories for each species. Only applicable if model="mixedRate" 
#' @param model The model of trait evolution (see details).
#' @param returnNodes Logical. If TRUE, alongside the tip values all node values are returned corresponding to APE's edge.matrix for the tree.
#' @param kappa Value of kappa transform.
#' @param lambda Value of lambda transform.
#' @param delta Value of delta transform.
#' @param alpha Value of alpha (OU) transform.
#' @param psi Value of psi transform.  Note that 'original nodes' from the full phylogeny can be included as an element on the phylogeny (e.g., phy$orig.node) as well as estimates of 'hidden' speciation (e.g., phy$hidden.speciation) if estimates of extinction (mu) are > 0.
#' @param lambda.sp Estimate of speciation (lambda) for the psi models
#' @param splitTime A split time (measured from the present, or most recent species) at which a shift in the rate occurs for the "timeSlice" model
#' @param timeRates The rates (from ancient to recent) for the timeSlice model
#' @param nodeIDs Integer - ancestral nodes of clades.
#' @param rateType If model="clade", a vector specifying if rate shift occurs in a clade ("clade") or on the single branch leading to a clade ("branch").
#' @param acdcRate Value of ACDC transform.
#' @param trend value of the expectation mean change through time
#' @param trend.anc.state the expected ancestal state for the trend model (default is 0)
#' @param branchLabels Branches on which different psi parameters are estimated in the "multipsi" model.
#' @param branchRates Numeric vector specifying relative rates for individual branches
#' @param rate a vector of relative rate parameters. The length of the vector is equal to the number of rates being estimated. 
#' @param group.means a vector of the relative difference in means between rate categories, expressed as a scalar applied to the expected standard deviation (see Ricklefs 2006)
#' @param cladeRates Numeric vector specifying telative rates for clades or logical to indicate scalar is included in the 'modeslice' model (the scalar is included in the mode.param argument with the 'modeslice' model).
#' @param rate.var Allows rate variation in BM modes in the 'modeslice' model
#' @param mode.order The order of modes for the 'modeslice' model. Any combination of 'BM', 'OU', 'acdc', and 'kappa'

#' @return Returns a matrix of simulated dated with taxon names as rownames (number of columns=n).
#' @references Ricklefs RE. 2006. Time, species, and the generation of trait variation in clades. Systematic Biology 55, 151-159.
#' @references Ricklefs RE. 2006. Thomas GH, Meiri S, & Phillimore AB. 2009. Body size diversification in Anolis: novel environments and island effects. Evolution 63, 2017-2030
#' @author Gavin Thomas, Mark Puttick
#' @seealso \code{\link{transformPhylo.ML}}, \code{\link{transformPhylo.ll}}, \code{\link{transformPhylo}}, \code{\link{transformPhylo.MCMC}}
#' @import mvtnorm
#' @examples
#' data(anolis.tree)
#' data(anolis.data)
#' 
#' # Simulate 10 sets of data with kappa=0.1 using the anolis tree
#' sim.dat1 <- transformPhylo.sim(phy=anolis.tree, n=10, model="kappa", kappa=0.1)
#' 
#' # Simulate 10 sets of data where rates and means differ between to the categories defined by "x"
#' x <- anolis.data$geo_ecomorph
#' names(x) <-  rownames(anolis.data)
#' sim.dat2 <- transformPhylo.sim(phy=anolis.tree, n=10, x=x, model="mixedRate", rate=c(1,1,2,4),
#' group.means=c(0,5,0,0))
#' @export
#' @import caper

transformPhylo.sim <- function(phy, n=1, x=NULL, model=NULL, returnNodes=FALSE, kappa=NULL, lambda=NULL, delta=NULL, alpha=NULL, psi=NULL, acdcRate=NULL, lambda.sp = NULL, trend=NULL, trend.anc.state=0, nodeIDs=NULL, rateType=NULL, cladeRates=NULL, branchRates=NULL, rate=NULL, group.means=NULL, splitTime=NULL, timeRates=NULL, branchLabels = NULL, rate.var=NULL, mode.order=NULL) {
	
	  model <- tolower(model)
	all.models <- c("bm", "trend", "kappa", "lambda", "delta", "free", "clade", "ou", "acdc", "psi", "multipsi", "timeslice", "mixedrate")
	if (any(is.na((match(model, all.models))))) stop(paste(model, "not recognised - please provide one of", paste0(all.models, collapse = ", ")))
	
	y.sim.function <- function(phy, return.node=FALSE) {	
		lengths <- phy$edge.length
		empty.mat <- matrix(NA, nrow=nrow(phy$edge), ncol=2)
		start.nodes <- which(phy$edge[,1] == Ntip(phy) + 1)
		empty.mat[start.nodes,1] <- rnorm(1)
	
		while(length(start.nodes) != 0) {
			next.node <- phy$edge[start.nodes,2]
			var.to.add <- rnorm(length(start.nodes), sd=sqrt(lengths[start.nodes]))
			empty.mat[start.nodes,2] <- empty.mat[start.nodes,1] + var.to.add
			start.nodes.int <- c()
			for(xx in 1:length(next.node)) {
				n.n <- which(phy$edge[,1] == next.node[xx])
				empty.mat[n.n,1] <- empty.mat[start.nodes[xx], 2]
				}
			start.nodes <- intersect(which(complete.cases(empty.mat[, 1])), which(!complete.cases(empty.mat[, 2])))
			}	
		tips <- which(phy$edge[, 2] <= Ntip(phy))
		y.out <- list()
		y.out$y <- matrix(empty.mat[tips, 2], dimnames=list(phy$tip.label))
		if(return.node) y.out$all.nodes <- empty.mat
		return(y.out)
		}

	switch(model,		  
		   
		   "bm" = {
					transformPhy <- phy
					ydum.all <- lapply(1:n, function(ix) y.sim.function(transformPhy, return.node=returnNodes))
					ydum <- sapply(ydum.all, function(re) re[[1]])
					rownames(ydum) <- phy$tip.label
					if(returnNodes) ydum <- list(trait.values=ydum, all.nodes=lapply(ydum.all, function(re) re[[2]]))
					},
			
			 "trend" = {
					transformPhy <- phy
					phyMat <- caper::VCV.array(transformPhy)
					attr(phyMat, "class") <- "matrix"
					tip.distance <- diag(vcv(transformPhy))
					trend.mean <- trend.anc.state + (tip.distance * trend)
					ydum <- as.matrix(t(mvtnorm::rmvnorm(n, mean=trend.mean, sigma = phyMat)))
					rownames(ydum) <- rownames(phyMat)
					if(returnNodes) warning("returnNodes not applicable to trend model, sorry")
					},
		   
		   "kappa" = {
					transformPhy <- transformPhylo(phy=phy, model="kappa", kappa=kappa, nodeIDs=nodeIDs)
					ydum.all <- lapply(1:n, function(ix) y.sim.function(transformPhy, return.node=returnNodes))
					ydum <- sapply(ydum.all, function(re) re[[1]])
					rownames(ydum) <- phy$tip.label
					if(returnNodes) ydum <- list(trait.values=ydum, all.nodes=lapply(ydum.all, function(re) re[[2]]))
					},
		   
		   "lambda" = {
					transformPhy <- transformPhylo(phy=phy, model="lambda", lambda=lambda)
					ydum.all <- lapply(1:n, function(ix) y.sim.function(transformPhy, return.node=returnNodes))
					ydum <- sapply(ydum.all, function(re) re[[1]])
					rownames(ydum) <- phy$tip.label
					if(returnNodes) ydum <- list(trait.values=ydum, all.nodes=lapply(ydum.all, function(re) re[[2]]))
					},
		   
		   "delta" = {
					transformPhy <- transformPhylo(phy=phy, model="delta", delta=delta, nodeIDs=nodeIDs)
					ydum.all <- lapply(1:n, function(ix) y.sim.function(transformPhy, return.node=returnNodes))
					ydum <- sapply(ydum.all, function(re) re[[1]])
					rownames(ydum) <- phy$tip.label
					if(returnNodes) ydum <- list(trait.values=ydum, all.nodes=lapply(ydum.all, function(re) re[[2]]))
					},
					
		   "free" = {
					transformPhy <- transformPhylo(phy=phy, model="free", branchRates=branchRates)
					ydum.all <- lapply(1:n, function(ix) y.sim.function(transformPhy, return.node=returnNodes))
					ydum <- sapply(ydum.all, function(re) re[[1]])
					rownames(ydum) <- phy$tip.label
					if(returnNodes) ydum <- list(trait.values=ydum, all.nodes=lapply(ydum.all, function(re) re[[2]]))
					},
		   
		   "clade" = {
					transformPhy <- transformPhylo(phy=phy, model="clade", nodeIDs=nodeIDs, cladeRates=cladeRates, rateType=rateType)
					ydum.all <- lapply(1:n, function(ix) y.sim.function(transformPhy, return.node=returnNodes))
					ydum <- sapply(ydum.all, function(re) re[[1]])
					rownames(ydum) <- phy$tip.label
					if(returnNodes) ydum <- list(trait.values=ydum, all.nodes=lapply(ydum.all, function(re) re[[2]]))
					},
		   
		   "ou" = {
		   	
		   			if (!is.ultrametric(phy)) {
      			        cophenetic.dist <- cophenetic.phylo(phy)
      			        vcv.matrix <- caper::VCV.array(phy)
      			        phyMat <- transformPhylo(phy=phy, model="OU", alpha=alpha, nodeIDs=nodeIDs, cophenetic.dist=cophenetic.dist, vcv.matrix=vcv.matrix)
      			        class(phyMat) <- "matrix"
      			        ydum <- as.matrix(t(mvtnorm::rmvnorm(n, sigma = phyMat)))
      			        rownames(ydum) <- rownames(phyMat)
      			        if(returnNodes) warning("returnNodes not applicable to OU model with non-ultrametric trees, sorry")
      			        } else {
      			        transformPhy <- transformPhylo(phy=phy, model="OU", alpha=alpha, nodeIDs=nodeIDs)
						ydum.all <- lapply(1:n, function(ix) y.sim.function(transformPhy, return.node=returnNodes))
						ydum <- sapply(ydum.all, function(re) re[[1]])
						rownames(ydum) <- phy$tip.label
						if(returnNodes) ydum <- list(trait.values=ydum, all.nodes=lapply(ydum.all, function(re) re[[2]]))
						}
					},
			"acdc" = {
					transformPhy <- transformPhylo(phy=phy, model="ACDC", acdcRate=acdcRate, nodeIDs=nodeIDs, cladeRates=cladeRates)
					ydum.all <- lapply(1:n, function(ix) y.sim.function(transformPhy, return.node=returnNodes))
					ydum <- sapply(ydum.all, function(re) re[[1]])
					rownames(ydum) <- phy$tip.label
					if(returnNodes) ydum <- list(trait.values=ydum, all.nodes=lapply(ydum.all, function(re) re[[2]]))
					},
		   "psi" = {
					transformPhy <- transformPhylo(phy = phy, model = "psi", psi = psi, lambda.sp = lambda.sp)
					ydum.all <- lapply(1:n, function(ix) y.sim.function(transformPhy, return.node=returnNodes))
					ydum <- sapply(ydum.all, function(re) re[[1]])
					rownames(ydum) <- phy$tip.label
					if(returnNodes) ydum <- list(trait.values=ydum, all.nodes=lapply(ydum.all, function(re) re[[2]]))
					},
			"multipsi" = {
       				transformPhy <- transformPhylo(phy = phy, model = "multipsi", psi = psi, lambda.sp = lambda.sp, branchLabels = branchLabels)
       				ydum.all <- lapply(1:n, function(ix) y.sim.function(transformPhy, return.node=returnNodes))
					ydum <- sapply(ydum.all, function(re) re[[1]])
					rownames(ydum) <- phy$tip.label
					if(returnNodes) ydum <- list(trait.values=ydum, all.nodes=lapply(ydum.all, function(re) re[[2]]))
    				},
					
			"timeslice" = {
				phy2 <- phy	   
		   		transformPhy <- transformPhylo(phy=phy, model="timeSlice", splitTime=splitTime, timeRates=timeRates)
		   		ydum.all <- lapply(1:n, function(ix) y.sim.function(transformPhy, return.node=returnNodes))
				ydum <- sapply(ydum.all, function(re) re[[1]])
				rownames(ydum) <- phy$tip.label
				if(returnNodes) ydum <- list(trait.values=ydum, all.nodes=lapply(ydum.all, function(re) re[[2]]))
		   		} ,
		   		
		   	"modeslice" = {
				phy2 <- phy	   
		   		transformPhy <- transformPhylo(phy=phy, model="modeSlice", mode.order=mode.order, splitTime=splitTime, cladeRates=cladeRates, rate.var=rate.var)
		   		ydum.all <- lapply(1:n, function(ix) y.sim.function(transformPhy, return.node=returnNodes))
				ydum <- sapply(ydum.all, function(re) re[[1]])
				rownames(ydum) <- phy$tip.label
				if(returnNodes) ydum <- list(trait.values=ydum, all.nodes=lapply(ydum.all, function(re) re[[2]]))
		   		} ,

			"mixedrate" = {
        			x <- as.matrix(x)
		        dat <- data.frame(x = x, y = rep(0, length(x[, 1])))
		        rateData <- as.rateData(y = "y", x = "x", rateMatrix = NULL, phy = phy, data = dat)
		        V <- transformRateMatrix(rateData, rate = rate)
		        expect.sd <- sqrt(mean(V[upper.tri(V)]))
		        if (is.null(group.means)) {
		            ydum.all <- lapply(1:n, function(ix) y.sim.function(transformPhy, return.node=returnNodes))
					ydum <- sapply(ydum.all, function(re) re[[1]])
					rownames(ydum) <- phy$tip.label
					if(returnNodes) ydum <- list(trait.values=ydum, all.nodes=lapply(ydum.all, function(re) re[[2]]))
		        } else {
		            x.means <- unique(rateData$x)
		            n.means <- length(x.means)
		            samp.means <- rep(NA, length(rateData$x))
		            ydum <- vector(mode = "list", length = length(group.means))
		            for (i in 1:n.means) {
		                samp.means[which(rateData$x == (i - 1))] <- rep(0 + (expect.sd * group.means[i]), length(which(rateData$x == (i - 1))))
		            }
		            ydum <- as.matrix(t(mvtnorm::rmvnorm(n, mean = samp.means, sigma = (V))))
		            rownames(ydum) <- rownames(V)
		       		if(returnNodes) warning("returnNodes not applicable to mixedrate model with non-ultrametric trees unique nodes, sorry")
		       		}
		       	}
		   )
	return(ydum)
}
