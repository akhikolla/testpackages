"read.indiv"<-function(file, skip=NA, cskip=NA){
  if(is.na(skip)){  skip <- getskip(file)}
  if(is.na(cskip)){cskip <- getcskip(file, skip)}
  IndivFileC <- scan(file, nlines=1, what="character",quiet=TRUE, skip=skip)
  if(cskip>0){IndivFileC <-IndivFileC[-(1:cskip)]}
  if(length(IndivFileC)==0){return(character(0))}
  if(length(IndivFileC)>=2 && IndivFileC[1]==IndivFileC[2]){
    IndivFileC <- IndivFileC[(1:length(IndivFileC))%%2==0]
  }  
  IndivFileC
}