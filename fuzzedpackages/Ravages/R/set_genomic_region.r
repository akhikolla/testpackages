set.genomic.region <- function(x, regions = genes.b37, flank.width = 0L) {
  if(!is.factor(regions$Gene_Name)) stop("regions$Gene_Name should be a factor with levels ordered in the genome order")

  # ce test est OK pour les facteurs aussi
  if(typeof(x@snps$chr) != "integer") 
    stop("x@snps$chr should be either a vector of integers, or a factor with same levels as regions$Chr")
  
  # remove duplicated regions if any
  w <- duplicated(regions$Gene_Name)
  if(any(w)) {
    regions <- regions[!w,]
  }

  # check if regions is sorted by chr / starting pos
  n <- nrow(regions)
  chr1 <- regions$Chr[1:(n-1)]
  chr2 <- regions$Chr[2:n]
  b <- (chr1 < chr2) | (chr1 == chr2 & regions$Start[1:(n-1)] <= regions$Start[2:n])
  if(!all(b)) {
    regions <- regions[ order(regions$Chr, regions$Start), ]
  }
  
  # if asked define larger regions
  if(is.finite(flank.width)) flank.width <- as.integer(flank.width)
  if(flank.width > 0L) {
    M <- as.integer(max(x@snps$pos, regions$End))  # joue le rôle de la position infinie !
    b <- regions$Chr[2:n] == regions$Chr[1:(n-1)]
    #Teste si les gènes se chevauchent
    b2 <- (regions$Start[2:n] > regions$End[1:(n-1)]) # non chevauchants
    start <- ifelse(b, ifelse(b2, as.integer(0.5*(regions$Start[-1] + regions$End[-n])), regions$Start[-1]),  0L)
    end <- ifelse(b, ifelse(b2, start-1L, regions$End[-n]), M)
    if(flank.width < Inf) {
      regions$Start <- pmax( c(0L, start), regions$Start - flank.width)
      regions$End <- pmin( c(end,M), regions$End + flank.width)
    } else {
      regions$Start <- c(0L,start)
      regions$End <- c(end,M)
    }
  }

  
  R <- .Call("label_multiple_genes", PACKAGE = "Ravages", regions$Chr, regions$Start, regions$End, x@snps$chr, x@snps$pos)
  R.genename <- unlist(lapply(R, function(z) paste(levels(regions$Gene_Name)[unlist(z)], collapse=",")))  
  R.genename[which(R.genename=="")] <- NA

  x@snps$genomic.region <- R.genename
  x@snps$genomic.region <- factor(x@snps$genomic.region, levels = unique(x@snps$genomic.region))
  x
}


