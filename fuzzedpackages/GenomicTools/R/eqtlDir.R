eqtlDir <- function(genoGroups,gex,mc=mc,nper, testType, MAF){

    innerFunction <- function(i)
    {
      # Calculate MAF
      maf.snp <- sum(genoGroups[genoGroups[,i]<3,i])/(2*(sum(genoGroups[,i]<3) ))
      maf.snp <- min(1-maf.snp, maf.snp)
      
      # Genotypes are coded as '3' if there is missing data, only proced if there is valid data
      if(sum(genoGroups[,i]==3)<(nrow(genoGroups))){
        missingData <- which((genoGroups[,i]==3)==TRUE)
      	PI <- length(table(genoGroups[genoGroups[,i]!=3,i]))
	    # All individuals have the same genotype, do nothing
	      if(PI==1){
	        innerOut <- 1.2
	
	      # Then the 2 groups comparison
	      } else if (PI==2){
	        if(testType=="asymptotic") testType <- "external"
	        innerOut <- gmw(gex[genoGroups[,i]!=3],genoGroups[genoGroups[,i]!=3,i],test="mw",type=testType,alternative="two.sided",nper=nper)$p.values

	      # And the three group comparison
	      } else if (PI==3){
	         p1 <- gmw(gex[genoGroups[,i]!=3],genoGroups[genoGroups[,i]!=3,i],test="triple",type=testType,alternative="greater",nper=nper)$p.values
	         p2 <- gmw(gex[genoGroups[,i]!=3],createGroups(genoGroups[genoGroups[,i]!=3,i],c(2,1,0)),test="triple",type=testType,alternative="greater",nper=nper)$p.values
	         innerOut <- min(2*min(p1,p2),1)
	  
	      }
      } else {
        innerOut <- 1.3
      }
      
      if(maf.snp<=MAF){
        innerOut <- 1.1
      }   
#      cat("InnerOut:", i, innerOut, "\n")     
      innerOut
    }
    
  if(is.matrix(genoGroups)){
    output <- unlist(mclapply(1:ncol(genoGroups),innerFunction,mc.cores=mc))    
  } else {
    genoGroups <- as.matrix(genoGroups)
    if(ncol(genoGroups)>1) genoGroups <- t(genoGroups)
    output <- unlist(innerFunction(1))
  }

  output
}
