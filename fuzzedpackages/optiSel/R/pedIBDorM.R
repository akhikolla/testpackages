
"pedIBDorM"<-function(Pedig, thisBreed=NA, keep.only=NULL, keep=keep.only){
  PedigAsDataTable <- "data.table" %in% class(Pedig)
  Pedig <- as.data.frame(Pedig)
  if(PedigAsDataTable){setDF(Pedig)}
  
  if(is.logical(keep)){keep<-Pedig[keep,1]}
  if(!is.null(keep)){keep<-as.character(keep); keep <- setdiff(keep, c(NA, "", " ", "0"))}

  Pedig <- prePed(Pedig, keep=keep)
  if(is.null(keep.only)){
    keep.only <- Pedig$Indiv
  }else{
    keep.only<-as.character(keep.only)
    keep.only<-Pedig$Indiv[Pedig$Indiv %in% keep.only]
  }
  
  
  Indiv<-1; Sire<-2; Dam<-3; Sex<-4; Breed<-5;
  if(is.na(thisBreed)){stop("The name of this breed is not specified.\n")}
  if(length(colnames(Pedig))==4){stop("Column breed is missing.\n")}
  
  for(i in c(Indiv, Sire, Dam, Breed)){Pedig[,i]<-as.character(Pedig[,i])}

  Rassen<-setdiff(names(table(Pedig[,Breed])),c(thisBreed))
  Selfing <- data.frame(
    paste(rep(c('Founder','Migrant'),20), rep(20:1,each=2),sep=''),
    c(NA, NA, paste(rep(c('Founder','Migrant'),19),rep(20:2,each=2),sep='')),
    c(NA, NA, paste(rep(c('Founder','Migrant'),19),rep(20:2,each=2),sep='')),
    NA, "Dummy", stringsAsFactors=FALSE)
  colnames(Selfing)<-colnames(Pedig)[1:5]
  
  Pedig<-rbind.data.frame(Selfing, Pedig[,1:5])

  Pedig[Pedig[,Breed] %in% Rassen, Sire] <- 'Migrant1'
  Pedig[Pedig[,Breed] %in% Rassen,  Dam] <- 'Migrant1'
  suppressWarnings(fOI <- 0.5*makeA(Pedig[,c(Indiv,Sire,Dam)])[keep.only, keep.only])
  #dimnames(fOI)<-list(Pedig[- (1:40),Indiv], Pedig[- (1:40),Indiv])
  Pedig[is.na(Pedig[,Sire])| Pedig[,Sire]=="0",Sire] <- 'Founder1'
  Pedig[is.na(Pedig[,Dam]) | Pedig[,Dam]=="0",  Dam] <- 'Founder1'
  Pedig[1:2, c(Sire, Dam)] <- NA
  #return(Pedig)  
  suppressWarnings(fII <- 0.5*makeA(Pedig[,c(Indiv,Sire,Dam)])[keep.only, keep.only]) #!!!
  #dimnames(fII)<-list(Pedig[- (1:40),Indiv], Pedig[- (1:40),Indiv])
  
  Res<-list()
  Res$pedIBDorM <- as(fOI + matrix(1,nrow=nrow(fII),ncol=ncol(fII)) - fII, "matrix")
  Res$pedIBDorMM<- as(fOI, "matrix")
  if(!is.null(keep.only)){for(i in names(Res))Res[[i]]<-Res[[i]][keep.only, keep.only]}
  class(Res)<-"kinMatrices"
  Res
}

