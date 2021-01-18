`print.vcf` <- function(x, n=6, m=6, fullHeader=FALSE,...){

  .Deprecated("GenomicTools.fileHandler::print.vcf", package="GenomicTools", msg="I/O Functions will be collected from now on in a new package GenomicTools.fileHandler")
  
# Correct for too large n and m
  nm <- dim(x$genotypes)
  nHeader <- min(n, length(x$header))
  n <- min(n,nm[1])
  m <- min(m,nm[2])

# Print the genotypes  
  cat("First",n,"rows and",m,"columns of $genotypes:\n")
  print(x$genotypes[1:n,x$map$snp.names[1:m], with=FALSE])
  cat("...",nrow(x$genotypes)-n,"rows and",ncol(x$genotypes)-m," columns omited \n\n")
  
# Print the header
  cat("First",n,"rows of $header:\n\n")  
  if(fullHeader){
    for(i in 1:nHeader){
      cat(x$header[i],"\n")
    }
  } else {
    for(i in 1:nHeader){
      cat(x$header[i],"\n")
    }
  }
  cat("...",length(x$header)-n,"rows omited \n\n")
  
  cat("First",n,"rows of $map:\n")  
  print(x$map[1:n,])
  cat("...",nrow(x$map)-n,"rows omited \n")
} 

