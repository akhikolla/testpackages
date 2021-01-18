"getcskip"<-function(file, skip){
  Line <- scan(file, nlines=1, what="character",quiet=TRUE, skip=skip+1)
  GT <- Line[length(Line)]
  cL   <- nchar(GT)
  if(cL==1){
    GT2<-Line[Line!=GT]
    GT2 <- GT2[length(GT2)]
    cskip <- length(Line)-which.max(rev(!c(FALSE,nchar(Line)==1 & Line %in% c(GT,GT2))))+1
  }else{
    if(cL!=3){stop("Genotype columns must contain 1 or 3 characters")}
    sep <- str_sub(GT, 2, 2)
    cskip <- length(Line)-which.max(rev(!c(FALSE,nchar(Line)==3 & str_sub(Line,2,2)==sep)))+1
  }
  cat(paste("Using cskip =",cskip,"\n"))
  cskip
}

