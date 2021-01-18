
"printParam4candes" <- function(cur, BreedNames){
  cur2      <- cur[cur$Type %in% c("trait", "kinship", "nat. kin."), c("Type", "Var",  "Breed", "obj.fun.and.constraints","Val")]
  cur2$Val  <- sprintf("%8.4f", cur2$Val)
  cur2$Breed[cur2$Breed %in% BreedNames] <- paste0("in ", cur2$Breed[cur2$Breed %in% BreedNames])
  cur2$Type <- paste0("for ", cur2$Type)
  cur2$Var  <- paste0(" '", cur2$Var, "' ")
  cur2$Type <- str_pad(as.character(cur2$Type),  max(nchar(cur2$Type)), "right")
  cur2$Var  <- str_pad(as.character(cur2$Var),   max(nchar(cur2$Var)), "right")
  cur2$Breed<- str_pad(as.character(cur2$Breed), max(nchar(cur2$Breed)), "right")
  if(identical(BreedNames, "missing")){cur2$Breed<-''}
  cur2$Type <- paste0(cur2$Type, cur2$Var, cur2$Breed, ":")
  cur2      <- cur2[,c("Type", "obj.fun.and.constraints", "Val")]
  cat("\n")
  colnames(cur2) <- c("Mean values of the parameters are:", "1", "  Value")
  print(cur2[,c(1,3),drop=FALSE], right = FALSE, row.names=FALSE)
  cat("\n")
  
  colnames(cur2) <- c("Available objective functions and constraints:", "2", "  Value")
  cur2[[1]] <- paste(cur2[[1]], cur2[[2]])
  print(cur2[,1,drop=FALSE], right = FALSE, row.names=FALSE, col.names=FALSE)
  cat("\n")
  cat(" ub  lb uniform\n\n")
}