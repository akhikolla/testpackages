
"pedIBDatN"<-function(Pedig, thisBreed=NA, keep.only=NULL, keep=keep.only, nGen=NA){
  getNe    <- !is.na(nGen)
  showInfo <- !is.null(keep.only)
  if("data.table" %in% class(Pedig)){
    Pedig <- as.data.frame(Pedig)
    setDF(Pedig)
  }
  if(is.na(thisBreed)){stop("The name of this breed is not specified.\n")}
  if(!("Indiv" %in% colnames(Pedig))){stop("Column Indiv is missing.\n")}
  if(!("Breed" %in% colnames(Pedig))){stop("Column Breed is missing.\n")}
  Pedig$Indiv <- as.character(Pedig$Indiv)
  Pedig$Breed <- as.character(Pedig$Breed)
  ids <- Pedig$Indiv[Pedig$Breed==thisBreed]
  
  if(is.logical(keep)){keep <- as.character(Pedig[keep,"Indiv"])}
  if(!is.null(keep)){  keep <- as.character(keep)}
  
  Pedig <- prePed(Pedig, keep=keep, lastNative=1234567, thisBreed=thisBreed)

  if(!is.null(keep)){
    keep <- Pedig$Indiv[Pedig$Indiv %in% keep]
    }
  
  if(is.null(keep.only)){
    keep.only <- ids
  }else{
    keep.only <- as.character(keep.only)
    keep.only <- ids[ids %in% keep.only]
  }
  
  if(("Sex" %in% colnames(Pedig)) && all(!is.na(Pedig[keep.only,"Sex"]))){
    sex <- Pedig[keep.only,"Sex"]
    if(all(sex=="male")){stop("Females are missing. Sexes can be ignored by removing column 'Sex'.\n")}
    if(all(sex=="female")){stop("Males are missing. Sexes can be ignored by removing column 'Sex'.\n")}
    if(!all(sex %in% c("female","male"))){stop("Sexes are not coded as 'female' and 'male'. They can be ignored by removing column 'Sex'.\n")}
  }
  
  if(getNe){
    Selection <- Pedig$Indiv
  }else{
    Selection <- keep.only
  }
  
  Rassen     <- setdiff(names(table(Pedig$Breed)), c(thisBreed))
  MigFounder <- Pedig$Indiv[(is.na(Pedig$Sire)|is.na(Pedig$Dam)) &   Pedig$Breed %in% Rassen]
  NatFounder <- Pedig$Indiv[(is.na(Pedig$Sire)|is.na(Pedig$Dam)) & !(Pedig$Breed %in% Rassen)]
  nMig <- length(MigFounder)
  nNat <- length(NatFounder)
  AMig <- matrix(2, nMig, nMig, dimnames=list(MigFounder, MigFounder))
  ANat <- matrix(2, nNat, nNat, dimnames=list(NatFounder, NatFounder))
  cat(paste0("Number of Migrant Founders: ", nrow(AMig), "\n"))
  cat(paste0("Number of Native  Founders: ", nrow(ANat), "\n"))
  cat(paste0("Individuals in Pedigree   : ", nrow(Pedig), "\n"))
  GB <- ((nrow(AMig)+nrow(ANat))^2 + length(Selection)^2 + nrow(Pedig)^2)*(7.45058066987776e-09)*1.1
  if((GB>1)&(length(Selection)>0.5*nrow(Pedig))){cat(paste0("Ensure that you have more than ", round(GB, 1), " GB memory available.\n"))}
  if(GB>1){cat("Computing fOI ...")}
  fOI  <- 0.5*makeA(Pedig[,1:3], AFounder=AMig, keep.only=Selection)[Selection, Selection]
  gc()
  if(GB>1){cat("finished\nComputing fII ...")}
  fII  <- 0.5*makeA(Pedig[,1:3], AFounder=adiag(AMig, ANat), keep.only=Selection)[Selection, Selection]
  if(GB>1){cat("finished\nCombining results ...")}
  rm(AMig)
  rm(ANat)
  gc()

  Cont       <- pedBreedComp(Pedig, thisBreed=thisBreed)
  NC         <- Cont[Selection, "native"]
  pedN       <- 1 - 0.5*(matrix(1 - NC, nrow=nrow(fII), ncol=ncol(fII), byrow=TRUE) +  matrix(1 - NC, nrow=nrow(fII), ncol=ncol(fII), byrow=FALSE)) - 0.5*(1-fII)
  pedIBDandN <- (fOI + 1 - fII) + pedN - 1
  dimnames(pedN)       <- list(Selection, Selection)
  dimnames(pedIBDandN) <- list(Selection, Selection)
  rm(fII)
  rm(fOI)
  pedN[pedN==0] <- 1e-14
  
  if(GB>1){cat("finished\n")}
  
  if(getNe){
    nativeNe   <- round(nativeNe(Pedig=Pedig, pedIBDandN=pedIBDandN, pedN=pedN, keep=keep.only, thisBreed=thisBreed, nGen=nGen),1)
    pedIBDandN <- pedIBDandN[keep.only, keep.only]
    pedN       <- pedN[keep.only, keep.only]
  }
  
  Res <- optiSolve::ratiofun(Q1=pedIBDandN, d1=0, Q2=pedN, d2=0, id=rownames(pedIBDandN))
  
  NC     <- Cont[rownames(pedIBDandN), "native"]
  Res$NC <- setNames(NC, rownames(pedIBDandN))
  
  if(showInfo){
    if(getNe){
      cat("Native effective size     : ", nativeNe, "\n", sep="")
      attr(Res,"nativeNe") <- nativeNe
    }
    Res$mean <- mean(pedIBDandN)/mean(pedN)
    #cat("Kinship at native alleles :",  round(Res$mean,     4), "\n")
    #cat("Native Contribution       : ", round(mean(Res$NC), 4), "\n", sep="")
  }else{
    Res$mean <- NA
  }
  
  return(Res)
}





nativeNe <- function(Pedig, pedIBDandN, pedN, keep, thisBreed, nGen=3){
  U  <- upper.tri(pedN[keep, keep])
  fD <- mean(pedIBDandN[keep, keep][U])/mean(pedN[keep, keep][U])
  ID <- keep
  
  for(i in 1:nGen){ID <- setdiff(unlist(Pedig[ID,c("Sire","Dam")]),c(NA,"0"))}
  ID <- intersect(Pedig$Indiv[Pedig$Breed==thisBreed],ID)
  U  <- upper.tri(pedN[ID, ID])
  fD0<- mean(pedIBDandN[ID, ID][U])/mean(pedN[ID, ID][U])
  1/(2*(1-((1-fD)/(1-fD0))^(1/nGen)))
}


