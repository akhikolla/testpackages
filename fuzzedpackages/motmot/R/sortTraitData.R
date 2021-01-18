#' @title Sort data and remove missing entries for tree and trait data
#' @description Plots a phylogeny with lines representing the value of a continuous trait
#' @param y A matrix of trait values with taxon names as rownames. Missing values should be NA
#' @param phy An object of class \code{phylo} or \code{multiPhylo}  (see \pkg{ape})
#' @param data.name If null the first column of y is assummed as the trait, otherwise if y is a matrix with more than one column either the name of the column or the number of the column must be supplied by data.name
#' @param log.trait Logical. If \code{TRUE}, data are log-transformed
#' @param pass.ultrametric Although trees that are believed to be ultrametric to pass the function \code{is.ultrametric} in \pkg{ape}
#' @return phy Tree with missing data pruned
#' @return trait Rearranged data with missing species removed
#' @author Mark Puttick
#' @examples 
#' data(anolis.tree)
#' data(anolis.data)
#' attach(anolis.data)
#' male.length <- matrix(Male_SVL, dimnames=list(rownames(anolis.data)))
#' any(is.na(male.length[,1]))
#' data.sorted <- sortTraitData(anolis.tree, male.length)
#' phy <- data.sorted[[1]]
#' male.length <- data.sorted[[2]]
#' @export

sortTraitData <- function (phy, y, data.name=NULL, log.trait = TRUE, pass.ultrametric=FALSE) 
{
	
	trait.names <- rownames(y)
	if(is.matrix(y) || is.data.frame(y)) {
		if(dim(y)[2] > 1) {
			if(dim(y)[2] > 1 && is.null(data.name)) stop("not sure which column to take - please provide a data.name")
			trait.data <- as.matrix(y[,data.name], dimnames=list(rownames(y)), nrow=nrow(y))
			} else {
			trait.data <- y	
			}
		} else {
			trait.data <- y
		}

    if (class(phy) == "multiPhylo") {
        tree <- phy
        phy <- phy[[1]]
        multi.phy <- TRUE
    } else {
        multi.phy <- FALSE
    }
    phy.names <- phy$tip.label
    missing <- apply(trait.data, 2, function(x) which(is.na(x)))
    missing <- c(unique(unlist(missing)))
    if(length(missing) > 0) {
    	dat.trait <- matrix(trait.data[-missing, ], ncol = ncol(trait.data), dimnames = list(trait.names[-missing]))
    } else {
    	dat.trait <- matrix(trait.data, ncol = ncol(trait.data), dimnames = list(trait.names))
    	}
    rm.tip <- match(phy.names, rownames(dat.trait))
    if (multi.phy) {
        red.phy <- lapply(tree, function(x) drop.tip(x, phy.names[which(is.na(rm.tip))]))
        class(red.phy) <- "multiPhylo"
        trait <- dat.trait[match(red.phy[[1]]$tip.label, rownames(dat.trait)), 
            ]
    } else {
        red.phy <- drop.tip(phy, phy.names[which(is.na(rm.tip))])
        trait <- dat.trait[match(red.phy$tip.label, rownames(dat.trait)), ]
    }
    
    if(pass.ultrametric == TRUE)  {
		outer <- red.phy$edge[,2]
		node.times <- nodeTimes(red.phy)
		externs <- which(outer <= Ntip(red.phy))
		red.phy$edge.length[externs] <- red.phy$edge.length[externs] + (node.times[externs,2] - 1e-7)
	}

    if (log.trait) trait <- log(trait)
    out <- list()
    out$phy <- red.phy
    out$trait <- as.matrix(trait)
    return(out)
}