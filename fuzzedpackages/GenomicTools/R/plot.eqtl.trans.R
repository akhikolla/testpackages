plot.trans <- function(x, genome){

  # Store the current par settings
   .pardefault <- par(no.readonly = T)
    par(mar=rep(0,4))
    circos.clear()
  
  # Adjust the input  
    genome$length <- genome$length/10^6
    gOrder <- as.character(genome$chr)
    genome$chr <- factor(genome$chr, levels=gOrder)

    
  # Basic circos graphic parameters
    circos.par(track.margin = c(0.1, 0))
    circos.par(cell.padding=c(0,0,0,0), track.margin=c(0,0.1), start.degree = 90, gap.degree =0)
    circos.par(points.overflow.warning=FALSE)
    
  # Chromosome details
    circos.initialize(factors = genome$chr, xlim = cbind(1, genome$length) )
  
  # Set the colors
    #cols.vals<-sapply(1:nrow(genome),function(x) c(runif(1),runif(1),runif(1)))
    rcols<-rainbow(nrow(genome))
    lcols<-rainbow(nrow(genome))
    
  # Plot chromosomes
    circos.trackPlotRegion(ylim = c(0, 1), factors = genome$chr, track.height=0.1,
                           # panel.fun for each sector
                             panel.fun = function(x, y) {
                             # Select details of current sector
                               name = get.cell.meta.data("sector.index")
                               i = get.cell.meta.data("sector.numeric.index")
                               xlim = get.cell.meta.data("xlim")
                               ylim = get.cell.meta.data("ylim")
                             
                               # Plot labels
                               circos.text(x=mean(xlim), y=1.8, labels=name, facing = "outside", cex=0.8)
                             
                               # Plot main sector
                                 circos.rect(xleft=xlim[1], ybottom=ylim[1], xright=xlim[2], ytop=ylim[2], col=rcols[i], border=rcols[i])
                             
                               # Blank in lower part of main sector
                                 circos.rect(xleft=xlim[1], ybottom=ylim[1], xright=xlim[2], ytop=ylim[1]+0.3, col = "white", border = "white")
                             
                               # White line all the way around
                                 circos.rect(xleft=xlim[1], ybottom=0.3, xright=xlim[2], ytop=0.32, col = "white", border = "white")
                             
                               # plot axis
                                # circos.axis(labels.cex=0.6, major.at=seq(from=0,to=floor(genome$length)[i],by=4), 
                                 #            labels.away.percentage = 0.15)
                               })
    

    wrongCHR <- FALSE
    linkedGenes <- as.character(unique(x$bed$Assoc.Gene))
    geneAnnot <- x$xAnnot[is.element(names(x$xAnnot),linkedGenes)]
    
    for(geneRun in 1:length(geneAnnot)){
      tmpAssoc <- x$bed[as.character(x$bed$Assoc.Gene)==names(geneAnnot)[geneRun],]
      tmpGene <- geneAnnot[[geneRun]]
      geneChr <- factor(tmpGene$Chr, levels=gOrder)
      
      for(assocRun in 1:nrow(tmpAssoc)){
        if(is.element(tmpAssoc$chr[assocRun],genome$chr)){
          snpChr <- factor(tmpAssoc$chr[assocRun], levels=gOrder)
          circos.link(sector.index1=geneChr, point1=c(tmpGene$Start/10^6,tmpGene$End/10^6),
                      sector.index2=snpChr, point2=c(tmpAssoc$Location[assocRun]/10^6),
                      col = lcols[geneChr])          
        } else {
          wrongCHR <- TRUE
        }

      }
    }
  if(wrongCHR) warning("There were sig. SNPs that coulnd't be matched to the genome information! Maybe wrong genome, or mislabeled chromosomes?")

# Restore the old par settings  
  par(.pardefault)  
}
