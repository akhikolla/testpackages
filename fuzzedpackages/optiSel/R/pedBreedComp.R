


"pedBreedComp"<-function(Pedig, thisBreed){
  PedigAsDataTable <- "data.table" %in% class(Pedig)
  Pedig <- as.data.frame(Pedig)
  if(PedigAsDataTable){setDF(Pedig)}
  ids <- Pedig$Indiv
  Pedig$Breed <- as.character(Pedig$Breed)

  Breeds <- setdiff(names(table(Pedig$Breed)), c(thisBreed))
  if("unknown"%in% Breeds){Breeds<-c(setdiff(Breeds, c("unknown")),"unknown")}
  SBreeds <- paste0("S.", Breeds)
  DBreeds <- paste0("D.", Breeds)
  
  Pedig[Pedig$Breed %in% Breeds & (is.na(Pedig[,2])|is.na(Pedig[,3])), 2] <- paste0("S.", Pedig$Breed[Pedig$Breed %in% Breeds  & (is.na(Pedig[,2])|is.na(Pedig[,3]))])
  Pedig[Pedig$Breed %in% Breeds & (is.na(Pedig[,2])|is.na(Pedig[,3])), 3] <- paste0("D.", Pedig$Breed[Pedig$Breed %in% Breeds  & (is.na(Pedig[,2])|is.na(Pedig[,3]))])

  Cont <- as.data.frame(genecont(Pedig[,1:3], from=c(SBreeds, DBreeds)))

  for(cName in setdiff(c(SBreeds, DBreeds), colnames(Cont))){
   Cont[,cName]<-0
  }
  Cont <- Cont[, SBreeds, drop=FALSE]+Cont[, DBreeds, drop=FALSE]
  colnames(Cont) <- Breeds
  Cont <- Cont[ids, rev(order(colMeans(Cont))), drop=FALSE]
  Cont <- data.frame(Indiv=rownames(Cont), native=1-rowSums(Cont), Cont, stringsAsFactors = FALSE)
  if(PedigAsDataTable){
    setDT(Cont)
    }
  Cont
}


