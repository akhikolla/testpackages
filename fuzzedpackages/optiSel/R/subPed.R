
"subPed"<-function(Pedig, keep, prevGen=3, succGen=0){
  PedigAsDataTable <- "data.table" %in% class(Pedig)
  Pedig <- as.data.frame(Pedig)
  if(PedigAsDataTable){setDF(Pedig)}
  colnames(Pedig)[1:3]<-c("Indiv", "Sire", "Dam")
  if(is.logical(keep)){keep<-Pedig$Indiv[keep]}
  Pedig<-prePed(Pedig, lastNative=1234567)
  selected  <- Pedig$Indiv %in% keep
  inPrevGen <- selected
  if(prevGen>0){
    for(i in 1:prevGen){
      inPrevGen <- inPrevGen | Pedig$Indiv %in% Pedig$Sire[inPrevGen] | Pedig$Indiv %in% Pedig$Dam[inPrevGen] 
    }
  }
  inSuccGen <- selected
  if(succGen>0){
    for(i in 1:succGen){
      inSuccGen <- inSuccGen |  Pedig$Sire %in% Pedig$Indiv[inSuccGen] |  Pedig$Dam %in% Pedig$Indiv[inSuccGen]
    }
    Sires <- Pedig$Sire[inSuccGen & !selected]
    Dams  <- Pedig$Dam[inSuccGen & !selected]
    inSuccGen <- inSuccGen | Pedig$Indiv %in% Sires | Pedig$Indiv %in% Dams 
  }  
  Pedig <- Pedig[selected | inPrevGen | inSuccGen, ]
  Pedig[!(Pedig$Sire %in% Pedig$Indiv), "Sire"] <- NA
  Pedig[!(Pedig$Dam %in% Pedig$Indiv),   "Dam"] <- NA
  Pedig<-prePed(Pedig, lastNative=1234567)
  Pedig$keep<-Pedig$Indiv %in% keep
  if(PedigAsDataTable){setDT(Pedig)}
  Pedig
}