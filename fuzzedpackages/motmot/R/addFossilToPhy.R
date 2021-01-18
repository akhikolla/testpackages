#' @title add a fossil to an interior branch of a time-scaled phylogeny
#'
#' @description the function takes a time-scaled phylogeny and adds a fossil at a branch in the past
#' @param phy Input phylogeny an object of class \code{phylo} (see \pkg{ape})
#' @param inGroup vector of tip labels defined as the ingroup - the fossil(s) will be placed on the stem branch leading to the 'inGroups' most recent common ancestor
#' @param fossil tip labels for the new fossil
#' @param fossilAge age of the fossil
#' @param minLength minimum length leading to the fossil 
#' @param maxLength maximum length leading to the fossil. If \code{NULL} (default) then the maximum bound is half the length of the branch leading to the crown node
#' @return the the time-scaled phylogeny with the fossil attached
#' @references Puttick, M. N., Kriwet, J., Wen, W., Hu, S., Thomas, G. H., & Benton, M. J. (2017). Body length of bony fishes was not a selective factor during the biggest mass extinction of all time. Palaeontology, 60, 727-741.
#' @author Mark Puttick
#' @examples
#' data(anolis.tree)
#' plot(anolis.tree)
#' nodelabels(214, 214)
#' # add fossil to node 214
#' in.groups <- node.descendents(x=214, phy=anolis.tree, tip.labels=TRUE)[[2]]
#' fossilPhy <- addFossilToPhy(anolis.tree, in.groups, fossil="fakeFossil", fossilAge=60)
#' plot(fossilPhy)
#' @export

addFossilToPhy <- function(phy, inGroup, fossil, fossilAge, minLength = 0.1, maxLength=NULL) {
	commonNode <- getMRCA(phy, inGroup)
	fossilTip <- read.tree(text="(fossilTip);")
	fossilTip$tip.label <- paste0(fossil)
	
	b.times <- nodeTimes(phy)
	edge.loc <- match(commonNode, phy$edge[,2])
	branchToCrown <- -diff(b.times[edge.loc,])
	minEdgeLength <- b.times[edge.loc, 2] -  fossilAge
	
	if(is.null(maxLength)) {
		maxLength <- branchToCrown
		} else {
		if(maxLength < branchToCrown) maxLength <- branchToCrown
	}
	
	if(minEdgeLength < 0) {
		minEdgeLength2 <- -minEdgeLength #+ minLength
		if(minEdgeLength2 > maxLength ) maxLength <- minEdgeLength2
		branchToCrown  <- runif(1, minEdgeLength2, maxLength)
		fossilTip$edge.length <- branchToCrown + b.times[edge.loc, 2] -  fossilAge
		return(bind.tree(phy, fossilTip, commonNode, branchToCrown))

	} else 
	{
		branchToCrown  <- runif(1, minLength, maxLength)
		fossilTip$edge.length <- branchToCrown + b.times[edge.loc, 2] -  fossilAge
		return(bind.tree(phy, fossilTip, commonNode, branchToCrown))
	}
}
