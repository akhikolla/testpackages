
makeA <- function(Pedig, keep.only=NULL, keep=keep.only, AFounder=NULL){
  PedigAsDataTable <- "data.table" %in% class(Pedig)
  if(PedigAsDataTable){
    Pedig <- as.data.frame(Pedig)
    setDF(Pedig)
    }
  
  if(is.logical(keep)){
    if(length(keep)!=nrow(Pedig)){stop("The length of logical vector keep must be equal to the number of rows in Pedig.\n")}
    keep <- Pedig[[1]][!is.na(keep)&keep] 
  }

  if(is.logical(keep.only)){
    if(length(keep.only)!=nrow(Pedig)){stop("The length of logical vector keep.only must be equal to the number of rows in Pedig.\n")}
    keep.only <- Pedig[[1]][!is.na(keep.only)&keep.only] 
  } 
  
  ids   <- as.character(Pedig[[1]])
  Pedig <- prePed(Pedig[,1:3], keep=keep, addNum=TRUE)

  if(is.null(keep.only)){
    keep.only <- ids
  }else{
    keep.only <- as.character(keep.only)
    keep.only <- setdiff(keep.only, c(NA, "", " ", "0"))
    keep.only <- ids[ids %in% keep.only]
  }

  if(is.null(AFounder)){
    AFounder   <- matrix(1, 0, 0)
    numFounder <- integer(0)
  }else{
    idF        <- Pedig$Indiv[Pedig$Indiv %in% rownames(AFounder)]
    AFounder   <- AFounder[idF, idF, drop=FALSE]
    numFounder <- Pedig[idF, "numIndiv"]
  }
  if(length(keep.only)<0.5*nrow(Pedig)){
    indKeep <- Pedig$Indiv[Pedig$Indiv %in% keep.only]
    numKeep <- Pedig[indKeep,"numIndiv"]
    Pedig$nOff <- 0
    x <- table(Pedig$Sire)
    Pedig[names(x),"nOff"] <- x
    x <- table(Pedig$Dam)
    Pedig[names(x),"nOff"] <- x
    pedKin  <- rcpp_makeA_lowMem(as.integer(Pedig$numSire), as.integer(Pedig$numDam), AFounder, as.integer(numFounder-1), as.character(indKeep), as.integer(numKeep-1), as.integer(Pedig$Indiv %in% indKeep), as.integer(Pedig$nOff))
  }else{
    pedKin <- rcpp_makeA(as.integer(Pedig$numSire), as.integer(Pedig$numDam), AFounder, as.integer(numFounder-1), as.character(Pedig$Indiv))
  }
  if(identical(Pedig$Indiv, keep.only)){return(pedKin)}
  pedKin[keep.only, keep.only]
}
