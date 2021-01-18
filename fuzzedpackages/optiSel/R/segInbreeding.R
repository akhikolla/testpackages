
"segInbreeding"<-function(files, map, minSNP=20, minL=1.0, unitP="Mb", unitL="Mb", a=0.0, keep=NULL, skip=NA, cskip=NA){
  ##################################################
  # Convert data tables to data frames             #
  ##################################################
  mapAsDataTable <- "data.table" %in% class(map)
  map <- as.data.frame(map)
  if(mapAsDataTable){setDF(map)}
  
  ##################################################
  #       Data preparation                         #
  ##################################################
  if(!("Name"  %in% colnames(map))){ stop("Column 'Name' is missing in map.")}
  if(!("Chr"   %in% colnames(map))){ stop("Column 'Chr'  is missing in map.")}  
  if(!is.null(keep)){
    keep <- as.character(keep)
    keep <- setdiff(keep, c(NA))
    }
  
  if(storage.mode(files)=="character"){
    files <- list(part1=files)
  }
  map$Chr <- as.character(map$Chr)
  for(i in 1:length(files)){
    names(files[[i]])<-str_sub(str_extract(basename(files[[i]]), "Chr[0-9,a-z]*"),4,-1)
  }
  if(!all(names(files[[1]]) %in% map$Chr)){
    stop("Some chromosomes are not in the map.")
  }
  if(any(duplicated(names(files[[1]])))){
    stop("For some chromosomes different files were provided.")
  }
  checkMap(map, unitL, unitP)
  rownames(map)<-map$Name
  if(is.na(skip)){  skip <- getskip(files[[1]][1])}
  if(is.na(cskip)){cskip <- getcskip(files[[1]][1], skip)}
  
  IndivFile1 <- scan(files[[1]][1], nlines=1, skip=skip, what="character", quiet=TRUE)
  if(cskip>0){IndivFile1 <-IndivFile1[-(1:cskip)]}
  if(length(IndivFile1)==1 || IndivFile1[1]!=IndivFile1[2]){
    IndivFile1 <- rep(IndivFile1, each=2)
  }
  if(is.null(keep)){
    index1 <- which(IndivFile1 %in% IndivFile1)
  }else{
    index1 <- which(IndivFile1 %in% as.character(keep))
    }
  NFile1 <- length(IndivFile1)
  Indiv1 <- IndivFile1[index1]
  N1     <- length(index1)
  if(N1==0){stop("No individuals are kept from file 1.")}
  
  if(length(files)>1){
    if(!identical(sort(names(files[[1]])), sort(names(files[[2]])))){
      stop("Chromosome names do not match for files 1 and files 2.\n")
    }
    IndivFile2 <- scan(files[[2]][1], nlines=1, skip=skip, what="character", quiet=TRUE)
    if(cskip>0){IndivFile2 <-IndivFile2[-(1:cskip)]}
    if(length(IndivFile2)==1 || IndivFile2[1]!=IndivFile2[2]){
      IndivFile2 <- rep(IndivFile2, each=2)
    }
    if(is.null(keep)){
      index2 <- which(IndivFile2 %in% IndivFile2)
    }else{
      index2 <- which(IndivFile2 %in% keep)
    }
    NFile2 <- length(IndivFile2)
    Indiv2 <- IndivFile2[index2]
    N2     <- length(index2)
  }else{
    IndivFile2 <- character(0)
    index2     <- numeric(0)
    NFile2     <- 0
    Indiv2     <- character(0)
    N2         <- 0
    files[[2]] <- rep("",length(files[[1]]))
  }  
  
  N <- N1 + N2
  Indiv <- c(Indiv1, Indiv2)
  symB <- scan(files[[1]][1], n=cskip+1, skip=1+skip, what="character", quiet=TRUE)[cskip+1]
  
  fROH  <- rep(0, N/2)
  gesL  <- 0
  
  for(chr in names(files[[1]])){
    cat(paste("Reading chromosome ", chr, "...  "))
    submap <- map[map$Chr==chr, ]
    M      <- nrow(submap)
    if(unitL %in% colnames(map)){cM <- submap[, unitL]}else{cM <- 1:M; cat("Using: unitL=Marker number\n")}
    if(unitP %in% colnames(map)){kb <- submap[, unitP]}else{kb <- 1:M; cat("Using: unitP=Marker number\n")}
    cM   <- (c(0,cM)+c(cM,cM[length(cM)]+cM[1]))/2
    kb   <- (c(0,kb)+c(kb,kb[length(kb)]+kb[1]))/2
    
    fROH <- fROH + rcpp_segInbreeding(as.character(files[[1]][chr]), as.character(files[[2]][chr]), as.integer(NFile1), as.integer(NFile2), as.integer(index1-1), as.integer(index2-1), as.integer(N1), as.integer(N2), as.integer(M), as.integer(minSNP), as.double(minL), as.double(cM), as.double(kb), as.double(a), as.character(symB), as.integer(skip), as.integer(cskip))
    gesL <- gesL + kb[length(kb)] - kb[1]
  }

  N   <- length(fROH)
  Res <- data.frame(Indiv=Indiv[2*(1:N)], Inbr=fROH/gesL, stringsAsFactors = FALSE)
  rownames(Res) <- Indiv[2*(1:N)]
  if(mapAsDataTable){setDT(Res)}
  Res
}
