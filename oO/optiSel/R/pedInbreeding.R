
"pedInbreeding"<-function(Pedig){
  PedigAsDataTable <- "data.table" %in% class(Pedig)
  Pedig <- as.data.frame(Pedig)
  if(PedigAsDataTable){setDF(Pedig)}
  
  Res <- data.frame(Indiv=Pedig[,1], Inbr=calcInbreeding(Pedig[,1:3]), stringsAsFactors = FALSE)
  rownames(Res) <- Pedig[,1]
  if(PedigAsDataTable){setDT(Res)}
  Res
}
  