#' @title Identify nodes and tips descended from a node 
#' @description Obtains a vector of the tips and nodes subtending from a node in a phylogeny.
#' @note \code{as.rateMatrix} calls the CAIC function \code{vcv.array} multiple times and this can be slow for large phylogenies (though faster than using the \pkg{ape} equivalent \code{vcv.phylo}).
#' @param x A positive integer
#' @param phy An object of class \code{phylo} (see \pkg{ape}).
#' @param tip.labels Logical - output tip.labels?
#' @details This function is stolen from the clade.members function in the CAIC package but returns both node and tip id's.
#' @return Returns a vector of node and tip ids descended from the tip(s) "x". If tip.labels=TRUE then returns a list of node ids and tip labels.
#' @author Gavin Thomas, David Orme
#' @export

node.descendents <- function (x, phy, tip.labels = FALSE) 
{
	 
	if (x <= Ntip(phy)) { return(numeric()) } else {
		node <- x
		Ntip <- Ntip(phy)
		ROOT <- Ntip + 1
		Nedge <- dim(phy$edge)[1]
		wbl <- !is.null(phy$edge.length)
		
		phy <- reorder(phy)
		root.node <- which(phy$edge[, 2] == node)
		start <- root.node + 1
		anc <- phy$edge[root.node, 1]
		next.anc <- which(phy$edge[-(1:start), 1] <= anc)
		keep <- if (length(next.anc)) start + 0:(next.anc[1] - 1) else start:Nedge
		
		nodes <- phy$edge[keep, 2]    
		nodes <- unique(nodes)
		if (tip.labels == TRUE) {
			nodesTips <- vector(mode = "list", length = 2)
			nodesTips[[1]] <- nodes[nodes > Ntip]
			nodesTips[[2]] <- with(phy, tip.label[nodes[nodes <=	Ntip]])
			return(nodesTips)
		}
		else {
			return(nodes)
		}}
}