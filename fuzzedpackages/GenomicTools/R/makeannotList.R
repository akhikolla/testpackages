makeAnnotList <- function(xAnnot){
  result <- list()
  if(sum(colnames(xAnnot)==c("Chr","Start","End","Gene"))!=4){
   stop("Please relabel the column names of xAnnot: 'Chr','Start','End', 'Gene'")
  }
  for(i in 1:nrow(xAnnot)){
    result[[i]] <- xAnnot[i,1:3] 
  }
  names(result) <- xAnnot[,4]
  result
}