
"gettraits"<-function(phen, quiet=FALSE){
  Trait       <- colnames(phen)
  singleBreed <- !(("Breed" %in% colnames(phen))&&(length(unique(phen$Breed))>1))
  colDesc     <- data.frame(" Trait"=Trait, validName=rep(TRUE,length(Trait)), isNumeric=TRUE, NA.ok=TRUE, row.names=seq_along(Trait), check.names = FALSE)

  for(i in seq_along(Trait)){
    if(Trait[i] %in% c("lb", "ub", "oc", "Born", "Sex", "Indiv", "Sire", "Dam", "Breed", "I", "Offspring", "herd", "Herd", "isCandidate", "c0","c1","Age","Class","Dead", "numIndiv", "numSire", "numDam")){
      colDesc[i,"validName"] <- FALSE
    }
    if(!is.numeric(phen[[i]])){
      colDesc[i,"isNumeric"] <- FALSE
    }
    if(all(is.na(phen[[i]]))){
      colDesc[i,"NA.ok"] <- FALSE
    }
    if(singleBreed && any(is.na(phen[[i]]))){
      colDesc[i,"NA.ok"] <- FALSE
    }
    if((!singleBreed) && any(is.na(phen[[i]]))){
      Breeds <- unique(phen$Breed[!is.na(phen[[i]])])
      if(length(Breeds)>1){
        colDesc[i,"NA.ok"] <- FALSE
      }
      if(any(is.na(phen[[i]]) & phen$Breed %in% Breeds)){
        colDesc[i,"NA.ok"] <- FALSE
      }
    }
  }
  
  use <- colDesc$validName & colDesc$isNumeric & colDesc$NA.ok
  
  if((!quiet) && any((!use) & colDesc$validName)){
    colDesc <- colDesc[colDesc$validName, ]
    cat("\n")
    cat("* Note:\n")
    print(colDesc,  row.names = FALSE)
    if(any(colDesc$validName & colDesc$isNumeric & !colDesc$NA.ok)){
      if(singleBreed){
        cat("\n")
        cat("  Traits must not contain NA values.\n")
        cat("\n")
      }else{
        cat("\n")
        cat("  Traits must contain values either \n")
        cat("  for exactly one breed or for all breeds.\n")
        cat("\n")
      }
    }
  }
  
  Trait[use]
}
