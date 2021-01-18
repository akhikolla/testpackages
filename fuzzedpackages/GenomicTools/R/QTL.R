# This script performs a QTL study for QTL studies 

QTL <- function(pheno, phenoSamples=NULL, geno=NULL, genoSamples=NULL, method="LM", mc=1, sig=NULL, testType="asymptotic", nper=2000, which=NULL, verbose=TRUE){

    # Initial values
      noNames <- FALSE
    
    # Test the case that the data is given in a one-row matrix (not casted as vector)
      if(is.matrix(pheno)){
        if(nrow(pheno)==1) pheno <- as.vector(pheno)
      }    
      
    # If only a vector of expression values is provided, transform it to a single row matrix
      if(is.null(pheno)) stop("No phenotype information provided in pheno!")
      if(is.null(geno)) stop("No genotypes provided in geno!")
      if(is.vector(pheno)){
        tmp <- names(pheno) 
        if(is.null(tmp)) noNames <- TRUE
        pheno <- as.matrix(pheno)
        rownames(pheno) <- tmp
      } else {
        if(is.null(rownames(pheno))) noNames <- TRUE
      }
    
    # Input checks
      method <- match.arg(method,c("LM","directional"))  
      testType <- match.arg(testType,c("permutation","asymptotic"))
    
    # Read in the genotype data if name is given, otherwise assume already imported SNPS are given as input
      if(is.character(geno)==TRUE && length(geno)==1)
      {
        # Check if there is a ped/map pair or vcf given as input
          fileEnding <- tolower(substr(geno, nchar(geno)-2,nchar(geno)))
          if(fileEnding=="ped"){
            case <- "ped"
            pedFile <- geno
            mapFile <- paste(substr(geno,1, nchar(pedFile)-3),"map",sep="")
          } else if(fileEnding=="vcf"){
            case <- "vcf"
            vcfFile <- paste(geno,".vcf",sep="")
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
      # CASE: No string provided, assume that genotype data was read in properly
       } else {
          if(class(geno)=="PedMap"){
             genoData <- geno
          } else if(class(geno)=="vcf"){
            genoData <- geno
            
      # In the other possibility is that the genotypes are just given in a matrix or vector form:    
          } else if(is.vector(geno)){
            # Here is just one vector with genotype information given
            if(length(geno)!=nrow(pheno)) stop("Amount of entered phenotypes and genotypes do not match!")
            genoData <- vectorToGenomatrix(geno)
      
          } else if(is.matrix(geno)||is.data.frame(geno)){ 
            # Here is a matrix with genotype information given
            if(nrow(geno)!=nrow(pheno)) stop("Amount of entered phenotypes and genotypes do not match!")
            genoData <- matrixToGenomatrix(geno)
            
          } else {
             stop("Please provide either a PedMap (importPED) of a VCF (importVCF) object, or the corresponding file path to either file.")
          }
       }

    
    # If no separate genoSamples object is given, we take those from the PedMap/VCF object:
      if(is.null(genoSamples)) genoSamples <-  rownames(genoData$genotypes)
    
    # Take only those genes from the annotation list that were requested
    if(!is.null(which)) {
      if(is.character(which)){
        pheno <- pheno[,is.element(colnames(pheno),which)]
      } else if(is.numeric(which)){
        oldSampleNames <- rownames(pheno)
        oldPhenoNames <- colnames(pheno)
        pheno <- pheno[,which]
      }
      if(length(which)==1){
        pheno <- matrix(pheno, ncol=1)
        colnames(pheno) <- oldPhenoNames[which]
        rownames(pheno) <- oldSampleNames
      }
    }
    
    if(noNames){
      if(is.null(phenoSamples)){
        if(verbose) cat("Phenotype are unnamed. We assume same order in phenotype and genotype objects and match samples based on that.\n")
        tempNames <- rownames(genoData$genotypes)
        rownames(pheno) <- tempNames
        phenoSamples <- tempNames
      } else {
        rownames(pheno) <- phenoSamples
      }
    }
    
    # In case that the row names have been changed, bring them into an order
      rownames(genoData$map) <- 1:nrow(genoData$map)
    
    # Sample statistics
      overlap <- is.element(rownames(pheno),genoSamples)
      olPerc <- round(sum(overlap)/nrow(pheno)*100,1)
      olPerc2 <- round(sum(overlap)/nrow(genoData$genotypes)*100,1)
      if(sum(overlap)==0) stop("No matching phenotype / genotype sample names!\n")
      if(verbose==TRUE) cat("We have for",olPerc,"% of the samples in the phenotype data the genotype information available. \n")
      if(verbose==TRUE) cat("We have for",olPerc2,"% of the samples in the genotype data the phenotype values available. \n")
    
    # Gene statistics
      if(is.null(colnames(pheno))) colnames(pheno) <- paste("Phenotype",1:ncol(pheno),sep="")
      matchingPheno <- colnames(pheno)
      if(verbose==TRUE) cat("We will test for",ncol(pheno),"phenotypes possible QTLs! \n")
      if(verbose==TRUE) cat("---------------------------------------------- \n")
      result <- list()
    
    # Reducing the phenotype data to those rows, where we have also genotype information available
      pheno <- pheno[is.element(rownames(pheno),genoSamples), ,drop=FALSE]
    
    # Now go through all possible genes
      qtl <- list()
    
    for(phenoRun in 1:length(matchingPheno)){
      qtlTemp <- list()
      SNPloc <- getSNPlocations(genotInfo=genoData$map, th=NULL)
        
      if(dim(SNPloc$SNPloc)[1]>0){
        SNPmatrix <- as.matrix(genoData$genotypes[,SNPloc$SNPloc$snp.names, with=FALSE])
        genoGroups <- matrix(as.numeric(SNPmatrix),nrow=nrow(SNPmatrix))
        colnames(genoGroups) <- colnames(SNPmatrix)
        rownames(genoGroups) <- rownames(genoData)
        genoGroups <- rearrange(genoGroups,rownames(pheno),genoSamples)
      
      #if(length(SNPloc$SNPloc$snp.names)>1){
      #  SNPmatrix <- as.matrix(genoData$genotypes[,SNPloc$SNPloc$snp.names])        
      #} else {
      #  SNPmatrix <- as.matrix(genoData$genotypes)
      #}

      #genoGroups <- matrix(as.numeric(SNPmatrix)-1,nrow=nrow(SNPmatrix))
      #colnames(genoGroups) <- colnames(SNPmatrix)
      #rownames(genoGroups) <- rownames(genoData)
      #genoGroups <- rearrange(genoGroups,rownames(pheno),genoSamples)
          
    # QTL case : LM
      if(method=="LM"){
        # if sig is set to Null all results will be reported - This might be very memory consuming!!!
          if(is.null(sig)){
              if(is.matrix(genoGroups)){
                qtlTemp[[phenoRun]] <- list(GeneLoc=rep(phenoRun,ncol(genoGroups)),TestedSNP=SNPloc[[1]],p.values=eqtlLM(genoGroups,pheno[,phenoRun], mc=mc))                
              } else {
                qtlTemp[[phenoRun]] <- list(GeneLoc=rep(phenoRun,1),TestedSNP=SNPloc[[1]],p.values=eqtlLM(genoGroups,pheno[,phenoRun], mc=mc))
              }

          } else {
              p.values <- eqtlLM(genoGroups,pheno[,phenoRun], mc=mc)
              pPos <- p.values<=sig
              qtlTemp[[phenoRun]] <- cbind(SNPloc[[1]][pPos,c(1,2,4)],p.values[pPos])
          }
        # QTL case: directional
        } else if(method=="directional"){
          # if sig is set to Null all results will be reported - This might be very memory consuming!!!
            if(is.null(sig)){ 
              if(is.matrix(genoGroups)){
              qtlTemp[[phenoRun]] <- list(GeneLoc=rep(phenoRun, ncol(genoGroups)),
                                          TestedSNP=SNPloc[[1]],
                                          p.values=eqtlDir(genoGroups,pheno[,phenoRun], mc=mc,nper=nper, testType=testType))
              } else {
                qtlTemp[[phenoRun]] <- list(GeneLoc=rep(phenoRun, 1),
                                            TestedSNP=SNPloc[[1]],
                                            p.values=eqtlDir(genoGroups,pheno[,phenoRun], mc=mc,nper=nper, testType=testType))
              }
            } else {
              p.values <- eqtlDir(genoGroups,pheno[,phenoRun],mc=mc,nper=nper, testType=testType)
              pPos <- p.values<=sig
              qtlTemp[[phenoRun]] <- cbind(SNPloc[[1]][pPos,c(1,2,4)],p.values[pPos])
            }
        }
      } else {
        warning("There were no variants within the window of gene ", matchingPheno[phenoRun])
      }
  
      
    # Join the output
      if(is.null(sig))
      {
        qtl[[phenoRun]] <- joinEQTL(qtlTemp)
      } else {
        if(sum(p.values<=sig)==0){
          qtl[[phenoRun]] <-  data.frame(chr="-1", SNP="-1", Location="-1", p.value="-1", Assoc.Pheno=matchingPheno[phenoRun], stringsAsFactors=FALSE)
        } else {
          #bedTemp <- joinEQTLsig(eqtlTemp)
          bedTemp <- do.call("rbind", qtlTemp)
          bedTemp <- cbind(bedTemp,rep(matchingPheno[phenoRun],max(1,nrow(bedTemp))))
          colnames(bedTemp) <- c("chr", "SNP", "Location", "p.value", "Assoc.Pheno")
          qtl[[phenoRun]] <- bedTemp          
        }
      }
      if(verbose==TRUE) cat ("We calculated QTLs for ",matchingPheno[phenoRun]," for ",prettyNum(nrow(SNPloc$SNPloc), big.mark = ",")," SNPs (",date(),")\n", sep="")
    }
    
    # Return the result
    if(is.null(sig))
    {
      names(qtl) <- matchingPheno
      result <- list(bed=NULL,qtl=qtl,pheno=pheno, geno=geno, phenoSamples=phenoSamples, genoSamples=genoSamples, method=method, mc=mc, which=which, type="full")
    } else {
      resBed <- do.call("rbind", qtl)
      remThese <- which(as.character(resBed[,1])== "-1")
      if(length(remThese)>0) resBed <- resBed[-remThese,]
      result <- list(bed= resBed, pheno=pheno, geno=geno, phenoSamples=phenoSamples, genoSamples=genoSamples, method=method, mc=mc, which=which, type="sig")
    }
    class(result) <- "qtlRes"
    result
}