genotypePlot <- function(snp, gene=NULL, eqtl=NULL, gex=NULL, geno=NULL, ylab=NULL, xlab=NULL, mainlab=TRUE){
  # Input checks
    if(is.null(eqtl) & is.null(gex) & is.null(geno)){
      stop("Please provide either an eqtl-object or expression and genotype objects")
    }
  # First bring the geno objects into the same order as the expression values are...
    if(is.null(eqtl)){
      eqtl$gex <- gex 
      eqtl$geno <- geno
    }
   if(!is.vector(eqtl$gex)){
     
   } else{
     if(sum(colnames(eqtl$gex)==gene)==0) stop("There are no expression values provided for the gene",gene)
   } 

   if(is.null(gene)) gene <- ""  

   rowGex <- rownames(eqtl$gex)
   takeThese <- is.element(as.character(eqtl$genoSamples),rowGex)
   eqtl$geno$genotypes <- eqtl$geno$genotypes[takeThese,]
   eqtl$geno$genotypes <- eqtl$geno$genotypes[order(rowGex),]
   snpCol <- which((eqtl$geno$map[,2]==snp)==TRUE)
   ..snpCol <- NULL  # This is just for the Cran check, hopefully it does not do any harm...
   snpValues <- as.numeric(as.vector(as.matrix(eqtl$geno$genotypes[, ..snpCol])))
   nGroups <- unique(snpValues)
   grExpr <- list()
   for(i in 1:length(nGroups)){
     #temp <- rownames(eqtl$geno$genotypes)[snpValues==i]
     temp <- eqtl$genoSamples[snpValues==i]
     if(gene==""){
       grExpr[[i]] <- eqtl$gex[is.element(rowGex,temp)]
     } else {
       grExpr[[i]] <- eqtl$gex[is.element(rowGex,temp),colnames(eqtl$gex)==gene]       
     }

   }
   
   temp <- eqtl$eqtl[names(eqtl$eqtl)==gene][[1]]
   temp2 <- temp[2]$TestedSNP[,2]==snp
   
   snpP <- temp[3]$p.values[temp2]
   
   ifelse(mainlab==TRUE, mainTitle <- paste(snp," and ",gene," (P-value:",snpP,")",sep=""), mainTitle <- "")
   
   boxplot(grExpr[[1]],xlim=c(0.5,length(grExpr)+0.5),ylim=c(min(as.vector(eqtl$gex)),max(as.vector(eqtl$gex))),main=mainTitle, ylab=ylab, xlab=xlab)
   if(length(grExpr)>1){
     for(i in 2:length(grExpr)){
       boxplot(grExpr[[i]],at=i,add=TRUE)
     }
   }  

   # Fill the x-axis
   tempRow <- eqtl$geno$map[eqtl$geno$map[,2]==snp,]
   al1 <- as.character(tempRow[5])
   al2 <- as.character(tempRow[6])
   axis(1,at=c(1,2,3),labels=c(paste(al1,"/",al1,sep=""),paste(al1,"/",al2,sep=""),paste(al2,"/",al2,sep="")))
}
