

"genecontALT" <- function(ID, Sire, Dam, NAncestors=NA) {
	ID    <- as.character(ID)
	Sire  <- as.character(Sire)
	Dam   <- as.character(Dam)
	Sire[is.na(Sire)] <- "0"
	Dam[is.na(Dam)]   <- "0"
  N     <- length(ID)
  if(is.na(NAncestors))NAncestors<-N
  
  Anteil<-matrix(0, nrow=N,ncol=NAncestors)
  colnames(Anteil)<-ID[1:NAncestors]
  rownames(Anteil)<-ID
  for(i in 1:N) {
    VaterAnt <-rep(0,NAncestors)
    MutterAnt<-rep(0,NAncestors)  
    if(sum(ID==Sire[i])>0)VaterAnt<-Anteil[Sire[i],]
    if(sum(ID==Dam[i])>0) MutterAnt<-Anteil[Dam[i],]
    Anteil[i,]<-(VaterAnt+MutterAnt)/2
    if(i<=NAncestors)Anteil[i,ID[i]]<-1
    }
  Anteil
}





