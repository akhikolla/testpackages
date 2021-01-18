#' @title prune tree and data to lineages present in a time bin in the past 
#'
#' @description the function takes a full tree and returns a pruned phylogeny with only tips and lineages found within a time bin preserved. If trait data are supplied the function will return tip states based either on the original tips found in the bin, or tip states inferred from ancestral states
#' @param phy An object of class \code{phylo} (see \pkg{ape}).
#' @param maxBin the start age (older time in myr from present) of the time bin in which lineages are preserved
#' @param minBin the final age (younger time in myr from present) of the time bin in which lineages are preserved
#' @param reScale if the most recent tip is not from the present, the age needed to add so the phylogeny is in 'real time'
#' @param allTraits a vector of trait data corresponding to the phy$edge object. The trait data represent tip and internal node data for the phylogeny
#' @param closest.min Logical. Should new tip values for lineages that span the bin be taken from the node nearest the \code{minBin} age (\code{closest.min=TRUE}, default) or the \code{maxBin} age \code{closest.min=FALSE}
#' @param traits.from.tip Logical. Should tip values for pendant edges in the bin be taken from the original tip value or the reconstructed node value (if it is closer than the tip value)
#' @return the pruned phylogeny. The object \code{descendants} refers to the lineages the branch in the time bin gave rise to before it was pruned. If traits are included a vector of trait values representing species at the tips.
#' @references Puttick, M. N., Kriwet, J., Wen, W., Hu, S., Thomas, G. H., & Benton, M. J. (2017). Body length of bony fishes was not a selective factor during the biggest mass extinction of all time. Palaeontology, 60, 727-741.
#' @author Mark Puttick
#' @examples
#' ## prune a random tree to taxa present between 4 and 2 units before present
#' # generate tree
#' set.seed(20)
#' tree <- rtree(20)
#' # generate traits
#' traits <- rnorm(20)
#' # plot tree and timeframe
#' plot(tree)
#' max.age <- nodeTimes(tree)[1,1]
#' abline(v=max.age - c(4, 2))
#' # prune tree to timeframe
#' cont.tree <- contemporaryPhy(phy=tree, maxBin=4, minBin=2, allTraits=traits)
#' plot(cont.tree$phy)
#' @export

contemporaryPhy <- function(phy, maxBin, minBin, reScale=0, allTraits, closest.min=TRUE, traits.from.tip=TRUE) {
	
	node.in <- match(unique(phy$edge[,1]), phy$edge[,1])
	node.times <- nodeTimes(phy)[node.in,1]
	names(node.times) <- unique(phy$edge[,1])
	
	startBin <- maxBin - reScale
	endBin <- minBin - reScale
	nodePreBin <- as.numeric(names(node.times)[which(node.times < endBin)])
	
	tip.times <- nodeTimes(phy)[which(phy$edge[,2] <= Ntip(phy)), 2]
	tipsInBin <- intersect(which(tip.times >= endBin), which(tip.times <= startBin))
	
	int.node <- which(phy$edge[,2] > Ntip(phy))
	descendantAsTrait <- rep(0, dim(phy$edge)[1])
	for(rr in int.node) descendantAsTrait[rr] <- length(node.descendents(phy$edge[rr,2], phy, T)[[2]])
	
	tipsInMat <- which(phy$edge[,2] <= Ntip(phy))
	terminalTips <- which(descendantAsTrait == 0)
	loseTips <- match(tipsInMat[tipsInBin], terminalTips)
	descendantAsTrait[terminalTips[-loseTips]] <- 1
	names(descendantAsTrait) <- phy$edge[,2]
	
	success <- timeTravelPhy(phy=phy, node=nodePreBin, nodeEstimate=descendantAsTrait, timeCut=endBin)
	suc <- removeNonBin(success$phy, success$tipData, keepByTime=startBin - endBin)
	
	out <- list()
	out <- suc
	names(out) <- c("phy", "descendants")
	
	if(!is.null(allTraits)) {
		
		first.times <- nodeTimes(phy)
		
		new.phylo <- suc$prunedPhy
		time.new.phylo <- nodeTimes(new.phylo) + reScale
	
		new.tip.names <- new.phylo$tip.label
		identify.tip <- regexpr("node_", new.phylo$tip.label)
		original.tip <- which(identify.tip == -1)
		same.as.original <- new.tip.names[original.tip]
		which.tips <- match(match(same.as.original, phy$tip.label), phy$edge[,2])
		new.spp <- which(identify.tip != -1)
		location.new.spp <- match(new.spp, new.phylo$edge[,2])
	
		node.original <- as.numeric(gsub("node_", "", new.phylo$tip.label[new.spp]))
		which.nodes <- match(node.original, phy$edge[,2])
		on.orig <- sort(c(which.nodes, which.tips))
	
		if (closest.min) closest.bin <- minBin else closest.bin <- maxBin
	
		trait.mat.new <- sapply(1:length(on.orig), function(i) {
			here.node <- on.orig[i]
			diff.age <- abs(first.times[here.node, ] - closest.bin)
			smallest.distance <- which.min(diff.age)
			if(smallest.distance == 2) {
				allTraits[here.node]
				} else {
				older.node <- match(phy$edge[here.node,1], phy$edge[,2])
				allTraits[older.node]
				}
			}
		)
		
		if(traits.from.tip)	trait.mat.new[original.tip] <- allTraits[which.tips]
		out$traits <- trait.mat.new
	}
	
	return(out)
}