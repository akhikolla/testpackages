
"conttac"<-function(cont, cohort, use=rep(TRUE,length(cohort)), mincont=0.05, long=TRUE){
  contAsDataTable <- "data.table" %in% class(cont)
  cont  <- cont[use==1,-1]
  cohort<- cohort[use==1]
  ContByYear <- stats::aggregate(cont,by=list(cohort),mean)
  rownames(ContByYear)<-ContByYear$Group.1
  native     <- ContByYear$native
  ContByYear <- ContByYear[,-(1:2)]
  Migrant    <- rowSums(ContByYear)
  ContByYear <- ContByYear[, colMeans(ContByYear)>=mincont, drop=FALSE]
  if(!("other" %in% colnames(ContByYear)))ContByYear$other<-0
  ContByYear$other <- ContByYear$other+Migrant-rowSums(ContByYear)
  if("unknown"%in% colnames(ContByYear)){
    ContByYear <- ContByYear[,c(setdiff(colnames(ContByYear), c("unknown")),"unknown")]
    }
  ContByYear <- cbind(native, ContByYear)
  ContByYear <- t(ContByYear)
  if(long){
    ContByYear <- melt(ContByYear)
    colnames(ContByYear) <- c("Breed","Year","Contribution")
    #ContByYear$Breed <- as.character(ContByYear$Breed)
  }
  if(contAsDataTable & long){setDT(ContByYear)}
  ContByYear
}

