eQTL <- function(gex=NULL, xAnnot=NULL, xSamples=NULL, geno=NULL, genoSamples=NULL, windowSize=0.5, method="directional", mc=1, sig=NULL, which=NULL, testType="asymptotic", nper=2000, verbose=TRUE, MAF=0.05 , IHaveSpace=FALSE){

  # Initial values
    noNames <- FALSE
    
  # Test the case that the data is given in a one-row matrix (not casted as vector)
    if(is.matrix(gex)){
      if(nrow(gex)==1) gex <- as.vector(gex)
    }  
  
  # Check if the xAnnot is in bed or in gtf format
    if(sum(class(xAnnot)=="gtf")>0){
      if(verbose) cat("xAnnot is given in gtf format (from importGTF). We transform it with gtfToBed() into required bed-format.\n")
      xAnnot <- gtfToBed(xAnnot)
    }
    
  # In case of a trans-eQTL set automatically a sig-value
    if(is.null(windowSize) & is.null(sig)){
      if(!IHaveSpace){
        sig <- 0.001
        warning("You choose trans-eQTL without specifying a 'sig'-value. This can lead to a large output, hence we set sig automatically to 0.001. If you want really
                 all results, please set the 'IHaveSpace'-Option to TRUE.")
      } else {
        if(verbose) cat("You set the IHaveSpace=TRUE option to get all results for the trans-eQTL run!\n")
      }
    } 
  
  # If only a vector of expression values is provided, transform it to a single row matrix
    if(is.null(gex)) stop("No gene expression provided in gex!")
    if(is.null(geno)) stop("No genotypes provided in geno!")
    if(is.vector(gex)){
      tmp <- names(gex) 
      oldNames <- names(gex)
      if(is.null(tmp)) noNames <- TRUE
      gex <- as.matrix(gex)
      if(verbose) cat("A vector of gene expression was provided. These expression values will be used for EACH gene in xAnnot. Please filter xAnnot accordingly, e.g. by using the 'which' option!\n")
      if(nrow(xAnnot)>1) gex <- gex[,rep(1,nrow(xAnnot))]
      colnames(gex) <- xAnnot[,4]
      rownames(gex) <- tmp
    } else {
      if(is.null(rownames(gex))) noNames <- TRUE
      oldNames <- rownames(gex)
    }

  # Input checks
    method <- match.arg(method,c("LM","directional"))  
    testType <- match.arg(testType,c("permutation","asymptotic"))
      
  # If the annotations are given as data frame, we will transform them into a list
    if(is.data.frame(xAnnot)){
      if(is.factor(xAnnot[,1])) xAnnot[,1] <- as.character(xAnnot[,1])
      if(is.factor(xAnnot[,4])) xAnnot[,4] <- as.character(xAnnot[,4])
      if(verbose==TRUE) cat("We will transform the gene annotations into a list ... ", sep="")
      xAnnot <- makeAnnotList(xAnnot)
      if(verbose==TRUE) cat("done (",date(),")!\n", sep="")
    }

  # Read in the genotype data if name is given, otherwise assume already imported SNPS are given as input
    if(is.character(geno)==TRUE)
    {
      # Check if there is a ped/map pair or vcf given as input
        fileEnding <- tolower(substr(geno, nchar(geno)-2,nchar(geno)))
        if(fileEnding=="ped"){
          case <- "ped"
          pedFile <- geno
          mapFile <- paste(substr(geno,1, nchar(pedFile)-3),"map",sep="")
        } else if(fileEnding=="vcf"){
          case <- "vcf"
          vcfFile <- geno
        } else {
          if(verbose) cat("No .ped or .vcf file ending detected in geno. Assume geno to be a ped/map filepair!\n")
          case <- "ped"
          pedFile <- paste(geno,".ped",sep="")
          mapFile <- paste(geno,".map",sep="")
        }
      # CASE: VCF
        if(case=="vcf"){
          if(verbose==TRUE) cat("Start reading the genotype information at",date(),"\n")
          if(verbose==TRUE) cat("vcf-file:",vcfFile,"\n")
          genoData <- importVCF(file=vcfFile)
      # CASE: PED/MAP
        } else if(case=="ped"){
          if(verbose==TRUE) cat("Start reading the genotype information at",date(),"\n")
          if(verbose==TRUE) cat("ped-file:",pedFile,"\n")
          if(verbose==TRUE) cat("map-file:",mapFile,"\n")
          genoData <- importPED(file=pedFile, snps=mapFile, verbose=FALSE)
        } 
    # Case: No string provided, assume that genotype data was read in properly
      } else {
        if(class(geno)=="PedMap"){
          genoData <- geno
        } else if(class(geno)=="vcf"){
          genoData <- geno
          
          # In the other possibility is that the genotypes are just given in a matrix or vector form:    
        } else if(is.vector(geno)){
          # Here is just one vector with genotype information given
       #   if(length(geno)!=nrow(pheno)) stop("Amount of entered phenotypes and genotypes do not match!")
          genoData <- vectorToGenomatrix(geno)
          
        } else if(is.matrix(geno)||is.data.frame(geno)){ 
          # Here is a matrix with genotype information given
      #    if(nrow(geno)!=nrow(pheno)) stop("Amount of entered phenotypes and genotypes do not match!")
          genoData <- matrixToGenomatrix(geno)
          
        } else {
          stop("Please provide either a PedMap (importPED) of a VCF (importVCF) object, or the corresponding file path to either file.")
        }
      }

  # Input checks
    if((ncol(gex)==1) & is.null(xAnnot))
    {
      warning("No annotations given, we will test all given SNPs against all given expressions!\n")
      xAnnot <- data.frame(Chr=as.character(genoData$map[,1]),Start=genoData$map[,4],End=genoData$map[,4],Gene=as.character(genoData$map[,2]))
      xAnnot <- makeAnnotList(xAnnot)
      windowSize=NULL
    }
    if((ncol(gex)>1) & is.null(xAnnot))
    {
      warning("No annotations given, we will test all given SNPs against all given expressions!\n gex is a matrix, hence this might take a while...\n")
      xAnnot <- data.frame(Chr=as.character(genoData$map[1,1]),Start=genoData$map[1,4],End=genoData$map[1,4], Gene="")
  #    xAnnot <- makeAnnotList(xAnnot)
      windowSize=NULL
    }

  # If no separate genoSamples object is given, we take those from the PedMap/VCF object:
    if(is.null(genoSamples)) genoSamples <-  rownames(genoData$genotypes)
  
  # Take only those genes from the annotation list that were requested
    if(!is.null(which)) {
      if(is.character(which)){
        xAnnot <- xAnnot[is.element(names(xAnnot),which)]
        gex <- gex[,is.element(colnames(gex),which)]          

      } else if(is.numeric(which)){
        xAnnot <- xAnnot[which]
        gex <- gex[,which]
      }
      
      if(length(which)==1){
        gex <- matrix(gex, ncol=1)
        colnames(gex) <- names(xAnnot)[1]
        rownames(gex) <- oldNames
      }
    }

    if(noNames){
     if(is.null(xSamples)){
       if(verbose) cat("Expression values are unnamed. We assume same order in expression and genotype objects and match samples based on that.\n")
       tempNames <- rownames(genoData$genotypes)
       rownames(gex) <- tempNames
       xSamples <- tempNames
     } else {
       rownames(gex) <- xSamples
     }
    }

  # In case that the row names have been changed, bring them into an order
    rownames(genoData$map) <- 1:nrow(genoData$map)

  # Sample statistics
    overlap <- is.element(rownames(gex),genoSamples)
    olPerc <- round(sum(overlap)/nrow(gex)*100,1)
    olPerc2 <- round(sum(overlap)/nrow(genoData$genotypes)*100,1)
    if(sum(overlap)==0) stop("No matching expression / genotype sample names!\n")
    if(verbose==TRUE) cat("We have for",olPerc,"% of the samples in the expression data the genotype information available. \n")
    if(verbose==TRUE) cat("We have for",olPerc2,"% of the samples in the genotype data the expression values available. \n")

  # Location statistics
    overlap <- is.element(colnames(gex),names(xAnnot))
    olPerc <- round(sum(overlap)/ncol(gex)*100,1)
    if(sum(overlap)==0) stop("No matching expression names / gene annotations!\n")
    if(verbose==TRUE) cat("We have for",olPerc,"% of the expression data the annotations. \n")

  # Gene statistics
    matchingGenes <- colnames(gex)[is.element(colnames(gex),names(xAnnot))]
    if(verbose==TRUE) cat("We will test for",length(matchingGenes),"genes possible eQTLs! \n")
    if(verbose==TRUE) cat("---------------------------------------------- \n")
    result <- list()

  # Reducing the expression data to those rows, where we have also genotype information available
    gex <- gex[is.element(rownames(gex),genoSamples), ,drop=FALSE]

  # Now go through all possible genes
    eqtl <- list()
 
    for(geneRun in 1:length(matchingGenes)){
    # Do that for each possible location of the gene (might not be unique...)
      tempAnnot <- xAnnot[[which((names(xAnnot)==matchingGenes[geneRun])==TRUE)]]
      eqtlTemp <- list()
   
      for(tempRun in 1:nrow(tempAnnot)){
      # SNP locations of variants inside the provided window
        SNPloc <- getSNPlocations(genotInfo=genoData$map, annot=tempAnnot[tempRun,], th=windowSize)
  
      # Run eQTL only if there are SNPs inside the window
        if(dim(SNPloc$SNPloc)[1]>0){
            SNPmatrix <- as.matrix(genoData$genotypes[,SNPloc$SNPloc$snp.names, with=FALSE])
            genoGroups <- matrix(as.numeric(SNPmatrix),nrow=nrow(SNPmatrix))
            colnames(genoGroups) <- colnames(SNPmatrix)
            rownames(genoGroups) <- rownames(genoData)
            genoGroups <- rearrange(genoGroups,rownames(gex),genoSamples)

          # eQTL case : LM
            if(method=="LM"){
            # if sig is set to Null all results will be reported - This might be very memory consuming!!!
    	        if(is.null(sig)){
    	            if(is.matrix(genoGroups)){
    	               eqtlTemp[[tempRun]] <- list(GeneLoc=rep(tempRun,ncol(genoGroups)),
    	                                           TestedSNP=SNPloc[[1]],
    	                                           p.values=eqtlLM(genoGroups,gex[,geneRun], mc=mc))
    	            } else {
    	              eqtlTemp[[tempRun]] <- list(GeneLoc=rep(tempRun,1),
    	                                          TestedSNP=SNPloc[[1]],
    	                                          p.values=eqtlLM(genoGroups,gex[,geneRun], mc=mc))
    	            }
    	        } else {
    	            p.values <- eqtlLM(genoGroups,gex[,geneRun], MAF=MAF, mc=mc)
              	  pPos <- p.values<=sig
    	            eqtlTemp[[tempRun]] <- cbind(SNPloc[[1]][pPos,c(1,2,4)],p.values[pPos])
    	        }
          # eQTL case: directional
            } else if(method=="directional"){
            # if sig is set to Null all results will be reported - This might be very memory consuming!!!
              if(is.null(sig)){ 
                 if(is.matrix(genoGroups)){
                   eqtlTemp[[tempRun]] <- list(GeneLoc=rep(tempRun, ncol(genoGroups)),
                                               TestedSNP=SNPloc[[1]],
                                               p.values=eqtlDir(genoGroups,gex[,geneRun], mc=mc, nper=nper, testType=testType, MAF=MAF))
                 } else {
                   eqtlTemp[[tempRun]] <- list(GeneLoc=rep(tempRun, 1),
                                               TestedSNP=SNPloc[[1]],
                                               p.values=eqtlDir(genoGroups,gex[,geneRun], mc=mc, nper=nper, testType=testType, MAF=MAF))
                   
                 }
  	        } else {
  	           p.values <- eqtlDir(genoGroups,gex[,geneRun],mc=mc,nper=nper, testType=testType, MAF=MAF)
  	           pPos <- p.values<=sig
  	           eqtlTemp[[tempRun]] <- cbind(SNPloc[[1]][pPos,c(1,2,4)],p.values[pPos])
          	}
          } 
       } else {
          warning("There were no variants within the window of gene ", matchingGenes[geneRun])
          p.values <- 2  # Take a p-value of 2 as internal indicator that there wasn't anything to test!
          eqtlTemp[[tempRun]] <- data.frame(GeneLoc=-1,TestedSNP=-1,p.values=-1, Assoc.Gene="-1")
       }
      }
      
      # Join the output
      if(is.null(sig))
      {
        eqtl[[geneRun]] <- joinEQTL(eqtlTemp)
        eqtl[[geneRun]]$GeneInfo <- tempAnnot
      } else {
        if(sum(p.values<=sig)==0){
          eqtl[[geneRun]] <-  data.frame(chr="-1", SNP="-1", Location="-1", p.value="-1", Assoc.Gene=matchingGenes[geneRun], stringsAsFactors=FALSE)
        } else {
          #bedTemp <- joinEQTLsig(eqtlTemp)
          bedTemp <- do.call("rbind", eqtlTemp)
          bedTemp <- cbind(bedTemp,rep(matchingGenes[geneRun],max(1,nrow(bedTemp))))
          colnames(bedTemp) <- c("chr", "SNP", "Location", "p.value", "Assoc.Gene")
          eqtl[[geneRun]] <- bedTemp          
        }
      }
    if(verbose==TRUE) cat ("We calculated eQTLs for ",matchingGenes[geneRun]," for ",prettyNum(nrow(SNPloc$SNPloc), big.mark = ",")," SNPs (",date(),")\n", sep="")
    }

  # Return the result
  if(is.null(sig))
  {
    names(eqtl) <- matchingGenes
    result <- list(bed=NULL,eqtl=eqtl,gex=gex, geno=geno, xAnnot=xAnnot, xSamples=xSamples, genoSamples=genoSamples, windowSize=windowSize, method=method, mc=mc, which=which, type="full")
  } else {
    resBed <- do.call("rbind", eqtl)
    remThese <- which(as.character(resBed[,1])== "-1")
    if(length(remThese)>0) resBed <- resBed[-remThese,]
    result <- list(bed= resBed,gex=gex, geno=geno, xAnnot=xAnnot, xSamples=xSamples, genoSamples=genoSamples, windowSize=windowSize, method=method, mc=mc, which=which, type="sig")
  }
  class(result) <- "eqtl"
  result
}



