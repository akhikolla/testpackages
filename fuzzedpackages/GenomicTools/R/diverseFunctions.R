# The QTL/eQTL function require an object that codes the genotypes as 00,01,02 (in $genotypes)
# and the locations in $map. If just a vector of genotypes is given, this function fakes the
# layout of it accoringly

vectorToGenomatrix <- function(x){
  out <- list(genotypes=NA,
              map=NA)
  x <- as.character(x)
  # The case that the genotypes are given in A/A notation
  if(grep("/",x[1])==1){
    alleles <- sort(unique(unlist(strsplit(as.character(unique(x)),"/"))))
    x.factor <- factor(as.character(x), levels = c(paste(alleles[1],"/",alleles[1],sep=""),
                                                   paste(alleles[1],"/",alleles[2],sep=""),
                                                   paste(alleles[2],"/",alleles[2],sep="")))
    levels(x.factor)[1] <- "01"
    levels(x.factor)[2] <- "02"
    levels(x.factor)[3] <- "03"
    out$genotypes <- as.matrix(x.factor)
    rownames(out$genotypes) <- 1:nrow(out$genotypes)
    colnames(out$genotypes) <- "inputSNP"
    
    out$map <- data.frame(snp.names="inputSNP",
                          allele.1=alleles[1],
                          allele.2=alleles[2])
    rownames(out$map) <- "inputSNP"
  }
  
  return(out)
}

matrixToGenomatrix <- function(x){
  sampleNames <- rownames(x)
  out <- list(genotypes=NA,
              map=NA)
  # The case that the genotypes are given in A/A notation
  for(i in 1:ncol(x)){
    if(grep("/",x[1,i])==1){
      alleles <- sort(unique(unlist(strsplit(as.character(unique(x[,i])),"/"))))
      x.factor <- factor(as.character(x[,i]), levels = c(paste(alleles[1],"/",alleles[1],sep=""),
                                                         paste(alleles[1],"/",alleles[2],sep=""),
                                                         paste(alleles[2],"/",alleles[1],sep=""),
                                                         paste(alleles[2],"/",alleles[2],sep="")))
      
      levels(x.factor)[levels(x.factor)==paste(alleles[1],"/",alleles[1],sep="")] <- "01"
      levels(x.factor)[levels(x.factor)==paste(alleles[1],"/",alleles[2],sep="")] <- "02"
      levels(x.factor)[levels(x.factor)==paste(alleles[2],"/",alleles[1],sep="")] <- "02"
      levels(x.factor)[levels(x.factor)==paste(alleles[2],"/",alleles[2],sep="")] <- "03"
      if(i==1){
        out$genotypes <- as.matrix(x.factor)        
      } else {
        tmp <- as.matrix(x.factor)
        out$genotypes <- cbind(out$genotypes, tmp)
      }
    }
    
    tmp <- data.frame(snp.names=colnames(x)[i],
                          allele.1=alleles[1],
                          allele.2=alleles[2])
    if(i==1){
      out$map <- tmp
    } else {
      out$map <- rbind(out$map,tmp)
    }
    rownames(out$map)[i] <- colnames(x)[i]
  }
  
  rownames(out$genotypes) <- sampleNames
  colnames(out$genotypes) <- colnames(x)
      
  return(out)
}

# Help function, taken from here: http://r.789695.n4.nabble.com/How-to-join-matrices-of-different-row-length-from-a-list-td3177212.html

cbind.fill <- function(...){
  nm <- list(...) 
  nm <- lapply(nm, as.matrix)
  n <- max(sapply(nm, nrow)) 
  do.call(cbind, lapply(nm, function (x) 
    rbind(x, matrix(, n-nrow(x), ncol(x))))) 
}