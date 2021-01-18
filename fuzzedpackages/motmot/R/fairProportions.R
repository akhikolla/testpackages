#' Calculate fair proportions phylogenetic diversity metric
#'
#' Calculate fair proportions phylogenetic diversity metric
#' Note that \code{as.rateMatrix} calls the CAIC function \code{vcv.array} multiple times and this can be slow for large phylogenies (though faster than using the "ape" equivalent \code{vcv.phylo}).
#' @param phy An object of class \code{phylo} (see \pkg{ape}).
#' @param nodeCount Logical - should root to tip node counts be returned (default is \code{FALSE})
#' @return Returns a matrix of fair proportion for all tips in phylogeny and node counts if selected.
#' @references Redding, D.W. and Mooers, A.O. (2006). Incorporating evolutionary measures into conservation prioritisation. Conservation Biology, 20, 1670-1678.
#' @references Isaac, N.J.B., Turvey, S.T., Collen, B., Waterman, C. and Baillie, J.E.M. (2007). Mammals on the EDGE: conservation priorities based on threat and phylogeny. PLoS ONE, 2, e296.
#' @author Gavin Thomas
#' @examples
#' data(anolis.tree)
#' 
#' fp <- fairProportions(anolis.tree)
#' fpNodes <- fairProportions(anolis.tree, nodeCount=TRUE)
#' @export
#' @import caper

fairProportions <- function (phy, nodeCount=FALSE) {
	
	treeMatrix <- caper::clade.matrix(phy)
	fpEdgeVals <- treeMatrix$edge.length / apply(treeMatrix$clade.matrix, 1, sum) 
	fpTips <- as.matrix(apply(fpEdgeVals * treeMatrix$clade.matrix, 2, sum))

	if (nodeCount == TRUE) {	
		nodeCount <- apply(treeMatrix$clade.matrix, 2, sum) - 1
		fpTips <- cbind(fpTips, nodeCount)
		colnames(fpTips) <- c("FP", "NodeCount")
		rownames(fpTips) <- phy$tip.label 
		} else {
		rownames(fpTips) <- phy$tip.label
		}

	rm(treeMatrix)
	return(fpTips)
}