#' @title Tree plotting for rates
#' @description Plots trees with colours based on rates of trait evolution. Also provides simple coloured plotting for trait values using the \code{\link[ape]{ace}} function in the \pkg{ape} library.
#' @param y A matrix of trait values.
#' @param x Output from \code{\link{summary.traitMedusa}}.
#' @param reconType Colour branches according to rate shifts ("rates" - requires traitMedusaObject) or ancestral state reconstruction ("picReconstruction"  - requires x).
#' @param palette Defines the colour scheme with four options: hotspot.colors (red to blue), heat.colors (yellow to red), cool.colors (blues), combi.colors (yellows to reds and blues)
#' @param ... Other functions to pass to \code{\link[ape]{plot.phylo}} 
#' @return Returns a data frame of colours used in plot along with rate (or ancestral state) range for each colour. 
#' @seealso \code{\link{transformPhylo.ML}}, \code{\link{summary.traitMedusa}}.
#' @author Gavin Thomas, Mark Puttick
#' @examples
#' # Data and phylogeny
#' data(anolis.tree)
#' data(anolis.data)
#'
#' # female SVL data
#' female.svl <- matrix(anolis.data[,"Female_SVL"], 
#' dimnames=list(rownames(anolis.data)))
#' input.data <- sortTraitData(phy=anolis.tree, y=female.svl, log.trait=TRUE)
#' 
#' # arbitarily reduce data size for speed in this example
#' phy.clade <- extract.clade(input.data[[1]], 182)
#' male.length.clade <- as.matrix(input.data[[2]][match(input.data[[1]]$tip.label, 
#' rownames(input.data[[2]])),])
#' # Identify rate shifts and print and plot results with up to one rate shifts 
#' # and minimum clade size of 10.
#' anolisSVL_MEDUSA <- transformPhylo.ML(male.length.clade, phy=phy.clade, 
#' model="tm1",minCladeSize=10, nSplits=1)
#' anolisSVL_MEDUSA_out <- summary(anolisSVL_MEDUSA, cutoff=1, AICc=FALSE)
#' colours <- plot(x = anolisSVL_MEDUSA_out,
#' reconType = "rates", type = "fan", cex=0.6, edge.width=3)
#' @import grDevices
#' @export plot.traitMedusa.model
#' @export

plot.traitMedusa.model <- function(x, y = NULL, ..., reconType="rates", palette="hotspot.colors") {
  
  if (class(x) != "traitMedusa.model")
	  stop("Please sapply output from transformPhylo.ML with model = 'tm1' or model = 'tm2'.")
	
  phy <- x$original.phy
  traitMedusaObject <- x
  
  cool.colors <- function (n, alpha = 1)
    grDevices::rainbow(n, start = 3/6, end = 4/6, alpha = alpha)
  heat.colors <- function (n, alpha = 1)
    grDevices::rainbow(n, start = 0, end = 1/6, alpha = alpha)
  hotspot.colors <- function (n, alpha = 1)
    grDevices::rainbow(n, start = 0, end = 4/6, alpha = alpha)
    
  switch (reconType,
    "picReconstruction" = {
       if (is.vector(y)==FALSE)
         stop("trait must be a vector with taxon names")
       a <- ace(y, phy, method="pic")
       reconTrait <- rep(NA, length(phy$edge[,2]))
       # internal nodes   
       reconTrait[match(names(a$ace)[2:length(a$ace)], phy$edge[,2])] <- a$ace[2:length(a$ace)]
       # tips
       idx <- match(names(y), phy$tip.label)
       reconTrait[match(idx, phy$edge[,2])] <- y
         
       ncolours <- Ntip(phy) + Nnode(phy) -1
       edgeColours <- rep(NA, ncolours)
         
       if (palette=="hotspot.colors")
         use.colors <- rev(hotspot.colors(ncolours))
       if (palette=="heat.colors")
         use.colors <- rev(heat.colors(ncolours))
       if (palette=="cool.colors")
         use.colors <- rev(cool.colors(ncolours))
       if (palette=="combi.colors")
         use.colors <- c(cool.colors(ncolours / 2), rev(heat.colors(ncolours / 2)))
         
       binsize <- (max(reconTrait) - min(reconTrait)) / length(reconTrait)
       color.bin <- matrix(NA, ncol=length(reconTrait), nrow=2)
       rownames(color.bin) <- c("min", "max")
         
       for (i in 1:length(reconTrait)) {
         color.bin[1, i] <- min(reconTrait) + (i - 1) * binsize
         color.bin[2, i] <- min(reconTrait) + i* binsize
         }
           
       color.bin[2,length(reconTrait)] <- max(reconTrait)	# necessary due to rounding problems
         
       for (i in 1: ncolours) {
         color.trait <- (reconTrait[i] >= color.bin["min",]) & (reconTrait[i] <= color.bin["max",])
         edgeColours[i] <- use.colors[which(color.trait==TRUE)]
       }
     },
		
	"rates" = {
	  if(traitMedusaObject$Rates[1,1] == "0")
	    stop("BM model - no shifts to plot")
	  cladeRates <- traitMedusaObject$Rates[,3]
	  nodes <- traitMedusaObject$Rates[,1]
	  rateType <- traitMedusaObject$Rates[,2]
	  ncolours <- Ntip(phy) + Nnode(phy) -1
	  nodeBranch <- data.frame(as.numeric(nodes), rateType)
	  nodes <- nodeBranch[nodeBranch[,"rateType"]=="clade",1]
	  branches <- nodeBranch[nodeBranch[,"rateType"]=="branch",1]
   	  cladeMembers <- cladeIdentity(phy, nodes)
	  edgeColours <- rep("black", length(phy$edge[,1]))
	
	  lgcladeRates <- log(as.numeric(cladeRates))
	  lgcladeRatesNodes <- lgcladeRates[nodeBranch[,"rateType"]=="clade"]
	  lgcladeRatesBranches <- lgcladeRates[nodeBranch[,"rateType"]=="branch"]
		
	  if (palette=="hotspot.colors")
	    use.colors <- rev(hotspot.colors(ncolours))
	  if (palette=="heat.colors")
	    use.colors <- rev(heat.colors(ncolours))
	  if (palette=="cool.colors")
	    use.colors <- rev(cool.colors(ncolours))
	  if (palette=="combi.colors")
	    use.colors <- c(cool.colors(ncolours / 2), rev(heat.colors(ncolours / 2)))
	  
	  binsize <- (max(abs(lgcladeRates)) - -max(abs(lgcladeRates))) / length(edgeColours)
	  color.bin <- matrix(NA, ncol=length(edgeColours), nrow=2)
	  rownames(color.bin) <- c("min", "max")
	  
	  for (i in 1:length(edgeColours)) {
	    color.bin[1,i] <- -max(abs(lgcladeRates)) + (i - 1) * binsize
	    color.bin[2,i] <- -max(abs(lgcladeRates)) + i * binsize
	    }
	  
	  color.bin[2, length(edgeColours)] <- max(abs(lgcladeRates))	# necessary due to rounding problems
	  color.bin[1, 1] <- -max(abs(lgcladeRates))	# necessary due to rounding problems

      for (i in 1:length(nodes)) {
      	color.trait <- (lgcladeRatesNodes[i] >= color.bin["min", ]) & (lgcladeRatesNodes[i] <= color.bin["max", ]) 
      	edgeColours[cladeMembers[,i]==1] <- use.colors[which(color.trait == TRUE)]
      	}
     
     for (i in 1:length(branches)) {
       color.trait <- (lgcladeRatesBranches[i] >= color.bin["min", ]) & (lgcladeRatesBranches[i] <= color.bin["max", ]) 
       edgeColours[which(phy$edge[,2]==branches[i])] <- use.colors[which(color.trait==TRUE)]
       }
     }
   )
   
   plot.phylo(phy, edge.color = edgeColours, ...)
   invisible(as.data.frame(t(rbind(use.colors, exp(color.bin)))))
   }