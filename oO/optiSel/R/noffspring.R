noffspring<-function(cand, N, random=TRUE){
  candAsDataTable <- "data.table" %in% class(cand)
  cand <- as.data.frame(cand)
  if(candAsDataTable){setDF(cand)}
  
  cand$Sex   <- as.character(cand$Sex)
  if(!("Indiv" %in% colnames(cand))){
    cand$Indiv <- rownames(cand)
  }
  cand$Indiv <- as.character(cand$Indiv)
  
  oc  <- cand$oc
  sex <- cand$Sex
  noff     <- rep(NA,length(oc))
  oc[oc<0] <- 0
  Sexes <- unique(sex)
  for(s in Sexes){
    noffs    <- floor(2*N*oc[sex==s])
    rest     <- 2*N*oc[sex==s]-noffs
    nMissing <- N-sum(noffs)
    if(nMissing>0){
      if(random){
        index  <- sample(1:sum(sex==s), size=nMissing, replace=FALSE, prob=rest)
      }else{
        index  <- rev(order(rest))[1:nMissing]
      }
      noffs[index] <- noffs[index]+1
    }
    noff[sex==s]<-noffs
  }
  
  Res <- data.frame(Indiv=cand$Indiv, nOff=noff, stringsAsFactors = FALSE)
  rownames(Res) <- cand$Indiv
  if(candAsDataTable){setDT(Res)}
  Res
}