eqtlLM.internal <- function(geno,gex, MAF){
  # Calculate MAF
    maf.snp <- sum(geno[geno<3])/(2*(sum(geno<3) ))
    maf.snp <- min(1-maf.snp, maf.snp)    
  
    res <- 1.3
    
    if(sum(is.na(geno)) < length(geno)){
      temp <- lm(gex~geno)
      ifelse(is.na(temp$coefficients[2]), res <- 1.2 ,  res <- summary(temp)$coefficients[2,4])
    } 
    
    if(maf.snp<=MAF) res <- 1.1
    
}

eqtlLM <- function(genoGroups, gex, MAF=MAF, mc=mc){
   if(is.matrix(genoGroups)){
      if(mc==1){
        res <- as.vector(apply(genoGroups,2,eqtlLM.internal, gex=gex, MAF=MAF))
      } else {
        res <- as.vector(unlist(mclapply(1:dim(genoGroups)[2], function(i) eqtlLM.internal(geno=genoGroups[,i], gex=gex, MAF=MAF), mc.cores=mc)))   
      }
    } else {
      genoGroups <- as.matrix(genoGroups)
      if(ncol(genoGroups)>1) genoGroups <- t(genoGroups)
      res <- as.vector(eqtlLM.internal(geno=genoGroups, gex=gex, MAF=MAF))
    }
  res
}

eqtlLM.naive <- function(genoGroups,gex){

  output <- c()
  for(i in 1:ncol(genoGroups))
  {
    if(sum(is.na(genoGroups[,i]))<(nrow(genoGroups))){
      fitData <- data.frame(x=genoGroups[,i],y=gex)
      temp <- lm(y~x,data=fitData)
      ifelse(is.na(temp$coefficients[2]), output[i] <- 1.25 ,  output[i] <- summary(temp)$coefficients[2,4])
    } else {
      output[i] <- 1.5
    }
  }
  output
}