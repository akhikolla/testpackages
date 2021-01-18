
"pedIBD"<-function(Pedig,  keep.only=NULL, keep=keep.only, kinFounder=NULL){
  AFounder <- NULL
  if(!is.null(kinFounder)){
    AFounder <- 2*kinFounder
  }
   0.5*makeA(Pedig, keep.only=keep.only, keep=keep, AFounder=AFounder)
}






