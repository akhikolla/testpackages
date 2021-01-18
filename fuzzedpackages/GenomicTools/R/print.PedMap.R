`print.PedMap` <- function(x, n=6, m=6, ...){

  .Deprecated("GenomicTools.fileHandler::printPedMap", package="GenomicTools", msg="I/O Functions will be collected from now on in a new package GenomicTools.fileHandler")
  
  # Correct for too large n and m
  nm <- dim(x$genotypes)
  n <- min(n,nm[1])
  m <- min(m,nm[2])
  
  cat("First",n,"rows and",m,"columns of $genotypes:\n")
  print(x$genotypes[1:n,x$map$snp.names[1:m], with=FALSE])
  cat("...",nrow(x$genotypes)-n,"rows and",ncol(x$genotypes)-m," columns omited \n\n")
  cat("First",n,"rows of $fam:\n\n")  
  print(x$fam[1:n,])
  cat("...",nrow(x$fam)-n,"rows omited \n\n")
  cat("First",n,"rows of $map:\n")  
  print(x$map[1:n,])
  cat("...",nrow(x$map)-n,"rows omited \n")
} 

