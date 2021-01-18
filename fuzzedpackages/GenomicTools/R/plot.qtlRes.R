`plot.qtlRes` <- function(x, which=NULL, sig=0.01, verbose=TRUE, log=FALSE, genome=NULL, pch=1, ...){
#plotQtlRes <- function(x, which=NULL, sig=0.01, verbose=TRUE, log=FALSE, genome=NULL, pch=1, ...){
  
  humanGenome68 <- data.frame(chr=c(1:22,"X","Y"),
                              length=c(249250621, 
                                       243199373, 
                                       198022430, 
                                       191154276, 
                                       180915260, 
                                       171115067, 
                                       159138663, 
                                       146364022, 
                                       141213431, 
                                       135534747, 
                                       135006516, 
                                       133851895, 
                                       115169878, 
                                       107349540, 
                                       102531392, 
                                       90354753, 
                                       81195210, 
                                       78077248, 
                                       59128983, 
                                       63025520, 
                                       48129895, 
                                       51304566, 
                                       155270560, 
                                       59373566)) 
  
   if(is.character(genome)){
     if(genome=="Human68") genome <- humanGenome68
   }

  if(is.null(genome)){
    warning("No genome information provided, we will visualize only the SNPs without further chromosomal length information!")  
    guessedChr <- unique(x$qtl[1][[1]]$TestedSNP$V1)
    gOrder <- as.character(guessedChr)
    nonNumeric <- which(is.element(tolower(gOrder),letters))
    if(length(nonNumeric)>0){
      chr <-as.numeric(gOrder[-nonNumeric])[order(as.numeric(gOrder[-nonNumeric]))]
      chr <- c(chr,sort(gOrder[nonNumeric]))      
    } else {
      chr <- as.numeric(gOrder)[order(as.numeric(gOrder))]
    }
    length <- numeric(length(chr))
    tmp <- x$qtl[1][[1]]$TestedSNP
    for(i in 1:length(chr)){
      length[i] <- max(tmp$V4[tmp$V1==chr[i]])
    }
    genome <- data.frame(chr=factor(chr, levels=chr),
                         length=length) 
    
  } else {
    # Adjust the input  
    genome$length <- genome$length
    gOrder <- as.character(genome$chr)
    genome$chr <- factor(genome$chr, levels=gOrder) 
  }


  if(is.null(which)) which <- 1:length(x$qtl)
  
# Visualize only significant loci  
  if(is.null(x$qtl)){

  } else {
# Visualize everything    
    for(phenoRun in which){
      tempQTL <- x$qtl[phenoRun][[1]]
      xLim <- c(1,sum(genome$length))
      if(log==TRUE){
        yLim <- c(0,-log(min(tempQTL$p.values)))
        yLab <- "-Log(p-values)"
      } else {
        yLim <- c(0,1)
        yLab <- "p-values"
      }

      offset <- 0
      plot(NA, xlim=xLim, ylim=yLim, ylab=yLab, xlab="Chromosomal Position", xaxt="n", main=colnames(x$pheno)[phenoRun])
      
      for(chrRun in 1:nrow(genome)){
        takeThese <- which(tempQTL$TestedSNP$V1==genome$chr[chrRun])
        tmpLocs <- tempQTL$TestedSNP$V4[takeThese]
        tmpP <- tempQTL$p.values[takeThese]
        if(log) tmpP <- -log(tmpP)
        tmpLocs <- tmpLocs + offset
        points(tmpLocs, tmpP, col=(chrRun%%2+1), pch=pch)
        offset <- offset + genome$length[chrRun]
      }

      axis(1, at=pairwiseDiffs(genome$length), labels=genome$chr)
      if(log) sig <- -log(sig)
      abline(h=sig, lty="dotted")
    } 
    
  }
}