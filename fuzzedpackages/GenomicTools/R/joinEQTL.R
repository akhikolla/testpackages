joinEQTL <- function(eqtlTemp){
  output <- list()
  # Join the ProbeLoc vectors
  temp <- c()
  for(i in 1:length(eqtlTemp))
  {
      temp <- c(temp,unlist(eqtlTemp[[i]][1]))
  }
  names(temp) <- c()
  output[[1]] <- temp

  # Join the testedGenome matrices
  #temp <- matrix(NA,ncol=6,nrow=0)
  temp <- setNames(data.table(matrix(nrow = 0, ncol = 6 )), colnames(eqtlTemp[[1]][2]$TestedSNP))
  for(i in 1:length(eqtlTemp))
  {
      if(is.element(-1,eqtlTemp[[i]][2]$TestedSNP)==TRUE){
        temp <- rbind(temp,  data.table(V1=-1,
                                        snp.names=-1,
                                        V3=-1,
                                        V4=-1,
                                        allele.1=-1,
                                        allele.2=-1))
      } else {
        temp <- rbind(temp,eqtlTemp[[i]][2]$TestedSNP)        
      }

  }
  output[[2]] <- as.data.frame(temp)


  # Join the p.values vectors
  temp <- c()
  for(i in 1:length(eqtlTemp))
  {
      temp <- c(temp,unlist(eqtlTemp[[i]][3]))
  }
  names(temp) <- c()
  output[[3]] <- temp
  names(output) <- c("ProbeLoc","TestedSNP","p.values")
  output

}
