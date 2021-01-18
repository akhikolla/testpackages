
"pedplot"<-function(Pedig, affected=NULL, status=NULL, label="Indiv", ...){
  PedigAsDataTable <- "data.table" %in% class(Pedig)
  Pedig <- as.data.frame(Pedig)
  if(PedigAsDataTable){setDF(Pedig)}
  id <- apply(Pedig[,label, drop=FALSE], 1, paste, collapse ="\n")
  if(is.null(affected) & ("keep" %in% colnames(Pedig))){affected<-Pedig$keep}
  if(is.null(status) & !is.null(affected) & ("Breed" %in% colnames(Pedig))){status<-!(Pedig$Breed%in%Pedig[affected>0,"Breed"])}
  Pedig$Sex[is.na(Pedig$Sex)]<-"unknown"
  Ped<-kinship2::pedigree(id=Pedig$Indiv, dadid=Pedig$Sire, momid=Pedig$Dam, sex=Pedig$Sex)
  kinship2::plot.pedigree(Ped, id=id, affected=affected, status=status, ...)
}